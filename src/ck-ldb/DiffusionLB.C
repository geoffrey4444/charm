/** \file DiffusionLB.C
 *  Authors: Monika G
 *           Kavitha C
 *
 */

/**
 *  1. Each node has a list of neighbors (bi-directional) (either topology-based
 *     or other mechanisms like k highest communicating nodes)
 *  2. Over multiple iterations, each node diffuses load to neighbor nodes
 *     by only passing load tokens (not actual objects)
 *  3. Once the diffusion iterations converge (load imbalance threshold is reached),
 *     actual load balancing is done by taking object communication into account
 */


#include "DiffusionLB.h"

#include "Diffusion.h"

#include "ck.h"
#include "ckgraph.h"
#include "envelope.h"
#include "elements.h"
#include "Heap_helper.C"
#define DEBUGF(x) CmiPrintf x;
#define DEBUGR(x) /*CmiPrintf x*/;
#define DEBUGL(x) CmiPrintf x;
#define NUM_NEIGHBORS 2
#define ITERATIONS 40

#define NODESIZE 1//2

// Percentage of error acceptable.
#define THRESHOLD 2

extern int getX(int node, int ny);
extern int getY(int node, int ny);

class Diffusion {
  public:
  void passPtrs(double *loadNbors,double *toSendLd,
                              double *toRecvLd, void (*func)(void*), void* obj);
  void setNeighbors(std::vector<int> neighbors, int neighborCount, double load);
};

/*readonly*/  CProxy_Diffusion diff_array;

//CreateLBFunc_Def(DiffusionLB, "The distributed graph refinement load balancer")
static void lbinit()
{
  LBRegisterBalancer<DiffusionLB>("DiffusionLB", "The distributed graph refine load balancer");
}

int CkMyNodeDiff() {
  int myPE = CkMyPe();
  return myPE/NODESIZE;
}

int CkNodeFirstDiff(int nd) {
  return nd*NODESIZE;
}

int CkNodeOfDiff(int myPE) {
  return myPE/NODESIZE;
}

int CkMyNodeSizeDiff() {
  return NODESIZE;
}

int CkNumNodesDiff() {
  return CkNumPes()/NODESIZE;
}

// TODO'S: 
// Topology
// Non migratable objects

void DiffusionLB::Strategy(const DistBaseLB::LDStats* const stats) {
  if (CkMyPe() == 0 && _lb_args.debug() >= 1) {
    double start_time = CmiWallTimer();
    DEBUGL(("In DiffusionLB strategy at %lf\n", start_time));
  }
  lb_started = false;
  AtSync();
}

DiffusionLB::DiffusionLB(CkMigrateMessage *m) : CBase_DiffusionLB(m) {
}

DiffusionLB::DiffusionLB(const CkLBOptions &opt) : CBase_DiffusionLB(opt) {
#if CMK_LBDB_ON
  lbname = "DiffusionLB";
  if (CkMyPe() == 0)
      DEBUGL(("[%d] Diffusion created\n",CkMyPe()));
  if (_lb_args.statsOn()) lbmgr->CollectStatsOn();
  //difflb_proxy = thisProxy;
  InitLB(opt);
#endif
}

DiffusionLB::~DiffusionLB()
{
#if CMK_LBDB_ON
  delete [] statsList;
  delete nodeStats;
  delete myStats;
  delete[] gain_val;
  delete[] obj_arr;
#endif
}

void DiffusionLB::InitLB(const CkLBOptions &opt) {
#if DEBUG_K
  DEBUGF(("[%d] InitLB\n", CkMyPe()));
#endif
  thisProxy = CProxy_DiffusionLB(thisgroup);
  if(CkMyPe()==0) {
    int NX = 2;
    int NY = CkNumNodesDiff()/NX;
    CProxy_BlockNodeMap map = CProxy_BlockNodeMap::ckNew(NX,NY, NODESIZE);
    CkArrayOptions opts(NX*NY);
    opts.setMap(map);
    diff_array = CProxy_Diffusion::ckNew(NX, NY, opts);
  }
  numNodes = CkNumNodesDiff();

  myspeed = lbmgr->ProcessorSpeed();
  // TODO: Initialize all class variables
  loadReceived = 0;
  statsReceived = 0;
  total_migrates = 0;
  total_migratesActual = -1;
  migrates_expected = -1;
  migrates_completed = 0;
  myStats = new DistBaseLB::LDStats;
  nodeFirst = CkNodeFirstDiff(CkMyNodeDiff());
  notif = 0;
#if DEBUG_K
  DEBUGL(("\n[PE-%d] nodeFirst = %d", CkMyPe(), nodeFirst));
#endif
  ComputeNeighbors();
  if(CkMyPe() == nodeFirst) {
    gain_val = NULL;
    obj_arr = NULL;
    my_load = 0;
    my_loadAfterTransfer = 0;
    statsList = new CLBStatsMsg*[nodeSize];
    nodeStats = new BaseLB::LDStats(nodeSize);
    loadPE.reserve(nodeSize);
    numObjects.reserve(nodeSize);
    //loadNeighbors = new double[neighborCount];//.reserve(neighborCount);
    prefixObjects.reserve(nodeSize);
    migratedTo.reserve(nodeSize);
    migratedFrom.reserve(nodeSize);
    for(int i = 0; i < nodeSize; i++) {
      loadPE[i] = 0;
      numObjects[i] = 0;
      migratedTo[i] = 0;
      migratedFrom[i] = 0;
      prefixObjects[i] = 0;
    }
  }
}

void DiffusionLB::AtSync() {
#if DEBUG_K
  DEBUGL(("\n[PE-%d]In DiffusionLB::AtSync()\n", CkMyPe())); fflush(stdout);
#endif
#if CMK_LBDB_ON
  if (!QueryBalanceNow(step()) || CkNumPes() == 1) {
    finalBalancing = 0;
    MigrationDone();
    return;
  }
  finalBalancing = 1;
  if(CkMyPe() == 0 && _lb_args.debug()) {
    receivedStats = 0;
  }
  migrates = 0;
  migratesNode = 0;
  objectHandles.clear();
  objectLoads.clear();
  // TODO: Check is it is the first load balancing step and then only 
  // perform this sending and QD
  if(step() == 0) {
    sendToNeighbors.clear();
    if(CkMyPe() == 0) {
      CkCallback cb(CkIndex_DiffusionLB::ProcessAtSync(), thisProxy);
      CkStartQD(cb);
    }
  }
  else {
    thisProxy[CkMyPe()].ProcessAtSync();
  }
#endif
}

void DiffusionLB::AddNeighbor(int node) {
#if 1//DEBUG_K
//  DEBUGL(("[PE-%d, Node-%d] My %luth neighbor node is %d\n", CkMyPe(), CkMyNodeDiff(), sendToNeighbors.size(), node));
#endif
  if(sendToNeighbors.size() > neighborCount)
    DEBUGL(("\n[PE-%d,node-%d]Adding nbors (node-%d) beyond count!! %lu>(of expected count %d)\n", CkMyPe(), CkMyNodeDiff(), node, sendToNeighbors.size(), neighborCount));
  sendToNeighbors.push_back(node);
}

void DiffusionLB::ProcessAtSync()
{
#if CMK_LBDB_ON
  start_lb_time = 0;
  internalBytes = externalBytes = received_nodes = 0;
#if DEBUG_K
  DEBUGL(("[%d] ProcessAtSync()", CkMyPe()));
#endif

  if (CkMyPe() == 0) {
    start_lb_time = CkWallTimer();
    if (_lb_args.debug())
      DEBUGL(("[%s] Load balancing step -%d starting at %f\n",
                              lbName(), step(), CkWallTimer()));
  }

#if DEBUG_K
  DEBUGL(("\nBefore AssembleStats()"));
#endif

  // assemble LB database
  statsmsg = AssembleStats();
  if(statsmsg == NULL)
    DEBUGL(("!!!Its null!!!\n"));

    // send to parent
//  DEBUGL(("[%d] Sending to parent #%d ReceiveStats\n", CkMyPe(), nodeFirst); fflush(stdout);
  CkMarshalledCLBStatsMessage marshmsg(statsmsg);
  thisProxy[nodeFirst].ReceiveStats(marshmsg);
  CkCallback cbm(CkReductionTarget(DiffusionLB, MaxLoad), thisProxy[0]);
  contribute(sizeof(double), &local_pe_load, CkReduction::max_double, cbm);
  CkCallback cba(CkReductionTarget(DiffusionLB, AvgLoad), thisProxy);
  contribute(sizeof(double), &local_pe_load, CkReduction::sum_double, cba);
  if(CkMyPe() != nodeFirst) {
    CkCallback cb(CkIndex_DiffusionLB::createNeighbors(), thisProxy);
    contribute(cb);
  }
}

void DiffusionLB::MaxLoad(double val) {
  DEBUGF(("\n[LB]Max PE load = %lf", val));
}

void DiffusionLB::AvgLoad(double val) {
  global_avg_load = val;
  if(CkMyPe()==0)
    DEBUGF(("\n[LB]Avg Node load = %lf", val/CkNumNodesDiff()));
}

void DiffusionLB::SumIntBytes(long val) {
  DEBUGF(("\nSum of Internal bytes = %lu", val));
}

void DiffusionLB::SumExtBytes(long val) {
  if(CkMyPe()==0)
    DEBUGF(("\nSum of External bytes = %lu", val));
//  thisProxy[CkMyPe()].ResumeClients(finalBalancing);
}

void DiffusionLB::doneNborExng() {
  if(CkMyPe() == CkNodeFirstDiff(CkMyNodeDiff()) && step() == 0) {
      //loadNeighbors.clear();
      //toReceiveLoad.clear();
      //toSendLoad.clear();
      sendToNeighbors.clear();
      loadNeighbors = new double[neighborCount];//.reserve(neighborCount);
      toReceiveLoad = new double[neighborCount];//.resize(neighborCount);
      toSendLoad = new double[neighborCount];//.resize(neighborCount);
      sendToNeighbors.reserve(neighborCount);
//      DEBUGL(("\n[PE-%d] setting sendToNeighbors size to %d", CkMyPe(), neighborCount);
      std::string nbor_nodes = " ";
      for(int i = 0; i < neighborCount; i++) {
        //Add your neighbors node-id as your neighbor
         nbor_nodes += "node-"+ std::to_string(nbors[i])+", ";
        AddNeighbor(nbors[i]);
      }
      sendToNeighbors.resize(neighborCount);
      DEBUGL(("\n[PE-%d,Node-%d] my neighbors: %s (#%d neighbors)\n", CkMyPe(), CkMyNodeDiff(), nbor_nodes.c_str(), neighborCount));
  }
  if(CkMyPe() == CkNodeFirstDiff(CkMyNodeDiff())) {
    int NX = 2;
    int NY = CkNumNodesDiff()/NX;
    Diffusion *diff_obj= diff_array(CkMyNodeDiff()).ckLocal();
    diff_obj->passPtrs(loadNeighbors, toSendLoad, toReceiveLoad, DiffusionLB::rLB, this);
    diff_obj->setNeighbors(sendToNeighbors, neighborCount,  my_load);
    
  //  thisProxy[CkMyPe()].iterate();
  }
#endif
}

void DiffusionLB::ComputeNeighbors() {
  //TODO: Use application graph for computing node neighbors
  //Assuming nodes are neighbors in a line
  nodeSize = CkMyNodeSizeDiff(); 
  nodeFirst = CkMyNodeDiff()*nodeSize;
  // TODO: Juan's topology aware mapping
  neighborCount = NUM_NEIGHBORS/2;
}

void DiffusionLB::sortArr(long arr[], int n, int *nbors)
{
 
  vector<std::pair<long, int> > vp;

  // Inserting element in pair vector
  // to keep track of previous indexes
  for (int i = 0; i < n; ++i) {
      vp.push_back(std::make_pair(arr[i], i));
  }

  // Sorting pair vector
  sort(vp.begin(), vp.end());
  reverse(vp.begin(), vp.end());

  int found = 0;
  for(int i=0;i<CkNumNodesDiff();i++)
    if(CkMyNodeDiff()!=vp[i].second) //Ideally we shouldn't need to check this
      nbors[found++] = vp[i].second;
  if(found == 0)
    DEBUGL(("\nPE-%d Error!!!!!", CkMyPe()));
}

// Assembling the stats for the PE
CLBStatsMsg* DiffusionLB::AssembleStats()
{
#if CMK_LBDB_ON
  // build and send stats
#if CMK_LB_CPUTIMER
  lbmgr->TotalTime(&myStats->total_walltime,&myStats->total_cputime);
  lbmgr->BackgroundLoad(&myStats->bg_walltime,&myStats->bg_cputime);
#else
  lbmgr->TotalTime(&myStats->total_walltime,&myStats->total_walltime);
  lbmgr->BackgroundLoad(&myStats->bg_walltime,&myStats->bg_walltime);
#endif
  lbmgr->IdleTime(&myStats->idletime);

  // TODO: myStats->move = QueryMigrateStep(step());

//    if (myStats->objData != NULL) { 
#if DEBUG_K
  DEBUGL(("Freeing \n"));
#endif
  myStats->objData.resize(lbmgr->GetObjDataSz());
  lbmgr->GetObjData(myStats->objData.data());

  myStats->commData.resize(lbmgr->GetCommDataSz());
  lbmgr->GetCommData(myStats->commData.data());

  const int osz = lbmgr->GetObjDataSz();
  const int csz = lbmgr->GetCommDataSz();

    // TODO: not deleted
  CLBStatsMsg* statsMsg = new CLBStatsMsg(osz, csz);
  statsMsg->from_pe = CkMyPe();

  // Get stats
#if CMK_LB_CPUTIMER
  lbmgr->GetTime(&statsMsg->total_walltime,&statsMsg->total_cputime,
                   &statsMsg->idletime, &statsMsg->bg_walltime,&statsMsg->bg_cputime);
#else
  lbmgr->GetTime(&statsMsg->total_walltime,&statsMsg->total_walltime,
                   &statsMsg->idletime, &statsMsg->bg_walltime,&statsMsg->bg_walltime);
#endif
  statsMsg->pe_speed = myStats->pe_speed;

  lbmgr->GetObjData(statsMsg->objData.data());
  lbmgr->GetCommData(statsMsg->commData.data());

  local_pe_load = 0;
  for (int i = 0; i < statsMsg->objData.size(); i++)
    local_pe_load += statsMsg->objData[i].wallTime;

  if(CkMyPe() == CkNodeFirstDiff(CkMyNodeDiff()))
    numObjects[0] = osz;
  return statsMsg;
#else
  return NULL;
#endif
}

void DiffusionLB::ReceiveStats(CkMarshalledCLBStatsMessage &&data)
{
#if CMK_LBDB_ON
  CLBStatsMsg *m = data.getMessage();
#if 1//DEBUG_K
//  DEBUGF(("[%d] GRD ReceiveStats from pe %d\n", CkMyPe(), m->from_pe));
#endif
  CmiAssert(CkMyPe() == nodeFirst);
  // store the message
  int fromRank = m->from_pe - nodeFirst;

  statsReceived++;
  AddToList(m, fromRank);

  if (statsReceived == nodeSize)  
  {
    // build LDStats
    BuildStats();
    if(step() == 0 && CkMyPe()==CkNodeFirstDiff(CkMyNodeDiff())) {
    //internalBytes = externalBytes = 0;
    long ebytes[CkNumNodesDiff()];
    std::fill_n(ebytes, CkNumNodesDiff(), 0);
    nbors = new int[NUM_NEIGHBORS+CkNumNodesDiff()];
    for(int i=0;i<CkNumNodesDiff();i++)
      nbors[i] = -1;
    neighborCount = NUM_NEIGHBORS/2;
//    DEBUGL(("\nedges = %lu", nodeStats->commData.size());
    for(int edge = 0; edge < nodeStats->commData.size(); edge++) {
      LDCommData &commData = nodeStats->commData[edge];
      if( (!commData.from_proc()) && (commData.recv_type()==LD_OBJ_MSG) ) {
        LDObjKey from = commData.sender;
        LDObjKey to = commData.receiver.get_destObj();

        // Check the possible values of lastKnown.
        int toPE = commData.receiver.lastKnown();
        int toNode = CkNodeOfDiff(toPE);
//        DEBUGL(("\ntoPE = %d",toPE);
        if(CkMyNodeDiff() != toNode && toNode!= -1) {
          ebytes[toNode] += commData.bytes;
      //    externalBytes += commData.bytes;
        } else if(CkMyNodeDiff() == toNode) {
      //    internalBytes += commData.bytes;
        }
      }
    }
//    for(int i=0;i<CkNumNodesDiff();i++)
//      DEBUGL(("\n[PE-%d,Node-%d] ebytes[to node %d] = %lu", CkMyPe(), CkMyNodeDiff(), i, ebytes[i]);
    sortArr(ebytes, CkNumNodesDiff(), nbors);
//    DEBUGL(("\n[PE-%d, node-%d], my largest comm neighbors are %d,%d\n", CkMyPe(), CkMyNodeDiff(), nbors[0], nbors[1]);
  }
  CkCallback cb(CkIndex_DiffusionLB::createNeighbors(), thisProxy);
  contribute(cb);
  }
}

void DiffusionLB::createNeighbors(){

    if(step()==0 && CkMyPe()==CkNodeFirstDiff(CkMyNodeDiff())) {
    for(int i=0;i<CkNumNodesDiff();i++) {
      int isNbor = 0;
      if(i != CkMyNodeDiff()){
        for(int j=0;j<NUM_NEIGHBORS/2;j++) {
          if(nbors[j] == i) {
            isNbor = 1;
            break;
          }
        }
      }
  #if 1//DEBUG
//      if(isNbor)
//        DEBUGL(("\n[PE-%d], notifying node %d [PE-%d]\n", CkMyPe(), i, CkNodeFirstDiff(i));
  #endif
      thisProxy[CkNodeFirstDiff(i)].notifyNeighbor(isNbor, CkMyNodeDiff());
    }
  } else
    doneNborExng();

    statsReceived = 0;

    // Graph Refinement: Generate neighbors, Send load to neighbors
//  }
#endif  
}

double DiffusionLB::average() {
  double sum = 0;
  DEBUGL(("\n[PE-%d load = %lf] n[0]=%lf, n[1]=%lf, ncount=%d\n", CkMyPe(), my_load, loadNeighbors[0], loadNeighbors[1], neighborCount));
  for(int i = 0; i < neighborCount; i++) {
    sum += loadNeighbors[i];
  }
  // TODO: check the value
  return (sum/neighborCount);
}
int DiffusionLB::GetPENumber(int& obj_id) {
  int i = 0;
  for(i = 0;i < nodeSize; i++) {
    if(obj_id < prefixObjects[i]) {
      int prevAgg = 0;
      if(i != 0)
          prevAgg = prefixObjects[i-1];
      obj_id = obj_id - prevAgg;
      break;
    }
  }
  return i;
}

bool DiffusionLB::AggregateToSend() {
  bool res = false;
  for(int i = 0; i < neighborCount; i++) {
    toSendLoad[i] -= toReceiveLoad[i];
    if(toSendLoad[i] > 0)
      res= true;
    toReceiveLoad[i] -= toSendLoad[i];
    DEBUGL(("[%d][node-%d] Diff: To Send load to node %d load %f res %d\n", CkMyPe(), CkMyNodeDiff(), sendToNeighbors[i], toSendLoad[i], res));
    DEBUGL(("[%d][node-%d] Diff: To Send load to node %d load %f res %d myLoadB %f\n", CkMyPe(), CkMyNodeDiff(), sendToNeighbors[i], toSendLoad[i], res, my_loadAfterTransfer));
  }
  return res;
}

void DiffusionLB::InitializeObjHeap(BaseLB::LDStats *stats, int* obj_arr,int n,
  int* gain_val) {
  for(int i = 0; i < n; i++) {
    obj_heap[i]=obj_arr[i];
    heap_pos[obj_arr[i]]=i;
  }
  heapify(obj_heap, ObjCompareOperator(&objs, gain_val), heap_pos);
}

int DiffusionLB::findNborIdx(int node) {
//  DEBUGL(("\n[PE-%d]Looking for node %d in nbor array of size %lu", CkMyPe(), node, sendToNeighbors.size());
  for(int i=0;i<sendToNeighbors.size();i++)
    if(sendToNeighbors[i] == node)
      return i;
  return -1;
}

#define SELF_IDX NUM_NEIGHBORS
#define EXT_IDX NUM_NEIGHBORS+1
void DiffusionLB::LoadBalancing() {
  int n_objs = nodeStats->objData.size();
  DEBUGL(("[%d] GRD: Load Balancing w objects size = %d \n", CkMyPe(), n_objs));

//  Iterate over the comm data and for each object, store its comm bytes
//  to other neighbor nodes and own node.

  //objectComms maintains the comm bytes for each object on this node
  //with the neighboring node
  //we also maintain comm within this node and comm bytes outside
  //(of this node and neighboring nodes)
  objectComms.reserve(n_objs);
  objectComms.resize(n_objs);

  if(gain_val != NULL)
      delete[] gain_val;
  gain_val = new int[n_objs];
  memset(gain_val, -1, n_objs);


  for(int i = 0; i < n_objs; i++) {
    objectComms[i].resize(NUM_NEIGHBORS+2);
    for(int j = 0; j < NUM_NEIGHBORS+2; j++)
      objectComms[i][j] = 0;
  }

  int obj = 0;
  for(int edge = 0; edge < nodeStats->commData.size(); edge++) {
    LDCommData &commData = nodeStats->commData[edge];
    // ensure that the message is not from a processor but from an object
    // and that the type is an object to object message
    if( (!commData.from_proc()) && (commData.recv_type()==LD_OBJ_MSG) ) {
      LDObjKey from = commData.sender;
      LDObjKey to = commData.receiver.get_destObj();
      int fromNode = CkMyNodeDiff();

      // Check the possible values of lastKnown.
      int toPE = commData.receiver.lastKnown();
      int toNode = CkNodeOfDiff(toPE);
      //store internal bytes in the last index pos ? -q
      if(fromNode == toNode) {
        int nborIdx = SELF_IDX;
        int fromObj = nodeStats->getHash(from);
        int toObj = nodeStats->getHash(to);
        //DEBUGR(("[%d] GRD Load Balancing from obj %d and to obj %d and total objects %d\n", CkMyPe(), fromObj, toObj, nodeStats->n_objs));
        objectComms[fromObj][nborIdx] += commData.bytes;
        // lastKnown PE value can be wrong.
        if(toObj != -1) {
          objectComms[toObj][nborIdx] += commData.bytes; 
        }
        internalBytes += commData.bytes;
      }
      else { // External communication
        int nborIdx = findNborIdx(toNode);
        if(nborIdx == -1)
          nborIdx = EXT_IDX;//Store in last index if it is external bytes going to non-immediate neighbors
        else {
          int fromObj = nodeStats->getHash(from);
          //DEBUGL(("[%d] GRD Load Balancing from obj %d and pos %d\n", CkMyPe(), fromObj, nborIdx);
          objectComms[fromObj][nborIdx] += commData.bytes;
          obj++;
        }
        externalBytes += commData.bytes;
      }
    }
    } // end for
    thisProxy[0].collectStats(CkMyNodeDiff(), internalBytes, externalBytes);
  }

void DiffusionLB::collectStats(int nodeId, long internalBytes_in, long externalBytes_in) {
  if(nodeId != CkMyNodeDiff()) {
    internalBytes += internalBytes_in;
    externalBytes += externalBytes_in;
  }
  received_nodes++;
  if(received_nodes == CkNumNodesDiff()) {
    CkPrintf("\n Total internalBytes = %lu, total externalBytes = %lu", internalBytes, externalBytes);
    for(int i=0;i<CkNumNodesDiff();i++)
      thisProxy[CkNodeFirstDiff(i)].continueLB();
  }
}

void DiffusionLB:: continueLB() {
    int n_objs = nodeStats->objData.size();

    // calculate the gain value, initialize the heap.
    double threshold = THRESHOLD*avgLoadNeighbor/100.0;

    actualSend = 0;
    balanced.resize(neighborCount);//toSendLoad.size());
    DEBUGL(("\nIterating through toSendLoad of size %lu", neighborCount));
    for(int i = 0; i < neighborCount; i++) {
      balanced[i] = false;
      if(toSendLoad[i] > 0) {
        balanced[i] = true;
        actualSend++;
      }
    }
    DEBUGL(("\n[PE-%d,Node-%d]actualSend = %d", CkMyPe(), CkMyNodeDiff(), actualSend));

    if(actualSend > 0) {

      if(obj_arr != NULL)
        delete[] obj_arr;

      obj_arr = new int[n_objs];

      for(int i = 0; i < n_objs; i++) {
        int sum_bytes = 0;
        //comm bytes with all neighbors
        vector<int> comm_w_nbors = objectComms[i];
        obj_arr[i] = i;
        //compute the sume of bytes of all comms for this obj
        for(int j = 0; j < comm_w_nbors.size(); j++)
            sum_bytes += comm_w_nbors[j];

        //This gives higher gain value to objects that have more within node communication
        gain_val[i] = 2*objectComms[i][SELF_IDX] - sum_bytes;
      }

      // T1: create a heap based on gain values, and its position also.
      obj_heap.clear();
      heap_pos.clear();
//      objs.clear();

      obj_heap.resize(n_objs);
      heap_pos.resize(n_objs);
 //     objs.resize(n_objs);
      std::vector<CkVertex> objs_cpy = objs;

      //Creating a minheap of objects based on gain value
      InitializeObjHeap(nodeStats, obj_arr, n_objs, gain_val); 

      // T2: Actual load balancingDecide which node it should go, based on object comm data structure. Let node be n
      int v_id;
      double totalSent = 0;
      int counter = 0;
      DEBUGL(("\n[PE-%d] my_loadAfterTransfer = %lf, actualSend=%d\n", CkMyPe(),my_loadAfterTransfer,actualSend));

      //return;
      while(my_loadAfterTransfer > 0 && actualSend > 0) {
        counter++;
        //pop the object id with the least gain (i.e least internal comm compared to ext comm)

        v_id = heap_pop(obj_heap, ObjCompareOperator(&objs, gain_val), heap_pos);
        
        DEBUGL(("\n On PE-%d, popped v_id = %d", CkMyPe(), v_id));
   
        /*If the heap becomes empty*/
        if(v_id==-1)          
            break;
        double currLoad = objs_cpy[v_id].getVertexLoad();
#if 1
        if(!objs[v_id].isMigratable()) {
          DEBUGL(("not migratable \n"));
          continue;
        }
#endif
        vector<int> comm = objectComms[v_id];
        int maxComm = 0;
        int maxi = -1;
        // TODO: Get the object vs communication cost ratio and work accordingly.
        for(int i = 0 ; i < neighborCount; i++) {
            
          // TODO: if not underloaded continue
          if(toSendLoad[i] > 0 && currLoad <= toSendLoad[i]+threshold) {
            if(i!=SELF_IDX && (maxi == -1 || maxComm < comm[i])) {
                maxi = i;
               maxComm = comm[i];
            }
          }
        }
        DEBUGL(("\n[PE-%d] maxi = %d", CkMyPe(), maxi));
          
        if(maxi != -1) {
#if 1
          migrates++;
          int pos = neighborPos[CkNodeOfDiff(nodeFirst)];
          int node = nbors[maxi];
          toSendLoad[maxi] -= currLoad;
          if(toSendLoad[maxi] < threshold && balanced[maxi] == true) {
            balanced[maxi] = false;
            actualSend--;
          }
          totalSent += currLoad;
          objs[v_id].setCurrPe(-1); 
          // object Id changes to relative position in PE when passed to function getPENumber.
          int objId = objs_cpy[v_id].getVertexId();
          if(objId != v_id) {
              DEBUGL(("\n%d!=%d", objId, v_id));fflush(stdout);
              CmiAbort("objectIds dont match \n");
          }
          int pe = GetPENumber(objId);
          migratedFrom[pe]++;
          int initPE = CkNodeFirstDiff(CkMyNodeDiff()) + pe;
          loadPE[pe] -= currLoad;
          numObjects[pe]--;
          DEBUGL(("[%d] GRD: Load Balancing object load %f to node %d and from pe %d and objID %d\n", CkMyPe(), currLoad, node, initPE, objId));
          // TODO: Change this to directly send the load to zeroth PE
          //thisProxy[nodes[node]].LoadTransfer(currLoad, initPE, objId);
          thisProxy[CkNodeFirstDiff(CkMyNodeDiff())].LoadMetaInfo(nodeStats->objData[v_id].handle, currLoad);
          thisProxy[initPE].LoadReceived(objId, CkNodeFirstDiff(node));
          my_loadAfterTransfer -= currLoad;
          int myPos = 0;//neighborPos[peNodes[nodeFirst]];
          loadNeighbors[myPos] -= currLoad;
          loadNeighbors[maxi] += currLoad;   
#endif
        }
        else {
          DEBUGL(("[%d] maxi is negative currLoad %f \n", CkMyPe(), currLoad));
        } 
      } //end of while
      DEBUGL(("[%d] GRD: Load Balancing total load sent during LoadBalancing %f actualSend %d myloadB %f v_id %d counter %d nobjs %lu \n",
          CkMyPe(), totalSent, actualSend, my_loadAfterTransfer, v_id, counter, nodeStats->objData.size()));
      for (int i = 0; i < neighborCount; i++) {
        DEBUGL(("[%d] GRD: Load Balancing total load sent during LoadBalancing toSendLoad %f node %d\n", CkMyPe(), toSendLoad[i], nbors[i]));
        }
      }//end of if
      // TODO: Put QD in intra node
      /* Start quiescence detection at PE 0.
      if (CkMyPe() == 0) {
          CkCallback cb(CkIndex_Diffusion::DoneNodeLB(), thisProxy);
          CkStartQD(cb);
      }*/
  DoneNodeLB();
  my_load = my_loadAfterTransfer;

}

// Load is sent from overloaded to underloaded nodes, now we should load balance the PE's within the node
void DiffusionLB::DoneNodeLB() {
  entered = false;
  if(CkMyPe() == nodeFirst) {
//    DEBUGL(("[%d] GRD: DoneNodeLB \n", CkMyPe()));
    double avgPE = averagePE();

    // Create a max heap and min heap for pe loads
    vector<double> objectSizes;
    vector<int> objectIds;
    minHeap minPes(nodeSize);
    double threshold = THRESHOLD*avgPE/100.0;
    
    for(int i = 0; i < nodeSize; i++) {
      if(loadPE[i] > avgPE + threshold) {
        DEBUGL(("[%d] GRD: DoneNodeLB rank %d is overloaded with load %f\n", CkMyPe(), i, loadPE[i]));
        double overLoad = loadPE[i] - avgPE;
        int start = 0;
        if(i != 0) {
          start = prefixObjects[i-1];
        }
        for(int j = start; j < prefixObjects[i]; j++) {
          if(objs[j].getCurrPe() != -1 && objs[j].getVertexLoad() <= overLoad) {
            objectSizes.push_back(objs[j].getVertexLoad());
            objectIds.push_back(j);
          }
        } 
      }
      else if(loadPE[i] < avgPE - threshold) {
        DEBUGL(("[%d] GRD: DoneNodeLB rank %d is underloaded with load %f\n", CkMyPe(), i, loadPE[i]));
        InfoRecord* itemMin = new InfoRecord;
        itemMin->load = loadPE[i];
        itemMin->Id = i;
        minPes.insert(itemMin);
      }
    }

    maxHeap objects(objectIds.size());
    for(int i = 0; i < objectIds.size(); i++) {
        InfoRecord* item = new InfoRecord;
        item->load = objectSizes[i];
        item->Id = objectIds[i];
        objects.insert(item); 
    }
    DEBUGR(("[%d] GRD DoneNodeLB: underloaded PE's %d objects which might shift %d \n", CkMyPe(), minPes.numElements(), objects.numElements()));

    InfoRecord* minPE = NULL;
    while(objects.numElements() > 0 && ((minPE == NULL && minPes.numElements() > 0) || minPE != NULL)) {
      InfoRecord* maxObj = objects.deleteMax();
      if(minPE == NULL)
          minPE = minPes.deleteMin();
      double diff = avgPE - minPE->load;
      int objId = maxObj->Id;
      int pe = GetPENumber(objId);
      if(maxObj->load > diff || loadPE[pe] < avgPE - threshold) {
          delete maxObj;
          continue;
      }
      migratedFrom[pe]++;
      DEBUGR(("[%d] GRD Intranode: Transfer obj %f from %d of load %f to %d of load %f avg %f and threshold %f \n", CkMyPe(), maxObj->load, pe, loadPE[pe], minPE->Id, minPE->load, avgPE, threshold));
      thisProxy[pe + nodeFirst].LoadReceived(objId, nodeFirst+minPE->Id);

      loadPE[minPE->Id] += maxObj->load;
      migratedTo[minPE->Id]++;
      loadPE[pe] -= maxObj->load;
      if(loadPE[minPE->Id] < avgPE) {
          minPE->load = loadPE[minPE->Id];
          minPes.insert(minPE);
      }
      else
          delete minPE;
      minPE = NULL;
    }
    while(minPes.numElements() > 0) {
      InfoRecord* minPE = minPes.deleteMin();
      delete minPE;
    }
    while(objects.numElements() > 0) {
      InfoRecord* maxObj = objects.deleteMax();
      delete maxObj;
    }

  // This QD is essential because, before the actual migration starts, load should be divided amongs intra node PE's.
    if (CkMyPe() == 0) {
      CkCallback cb(CkIndex_DiffusionLB::MigrationEnded(), thisProxy);
      CkStartQD(cb);
    }
    /*for(int i = 0; i < nodeSize; i++) {
      thisProxy[nodeFirst + i].MigrationInfo(migratedTo[i], migratedFrom[i]);
    }*/
  }
}

double DiffusionLB::averagePE() {
  int size = nodeSize;
  double sum = 0;
  for(int i = 0; i < size; i++) {
      sum += loadPE[i];
  }
  return (sum/(size*1.0)); 
}

int DiffusionLB::FindObjectHandle(LDObjHandle h) {
  for(int i = 0; i < objectHandles.size(); i++)
    if(objectHandles[i].id == h.id)
      return i;
  return -1;  
}

void DiffusionLB::LoadReceived(int objId, int fromPE) {
  // load is received, hence create a migrate message for the object with id objId.
  if(objId > myStats->objData.size()) {
    DEBUGR(("[%d] GRD: objId %d total objects %d \n", objId, myStats->objData.size()));
    CmiAbort("this object does not exist \n");
  }
  MigrateInfo* migrateMe = new MigrateInfo;
  migrateMe->obj = myStats->objData[objId].handle;
  migrateMe->from_pe = CkMyPe();
  migrateMe->to_pe = fromPE;
  //migrateMe->async_arrival = myStats->objData[objId].asyncArrival;
  migrateInfo.push_back(migrateMe);
  total_migrates++;
  entered = false;
  DEBUGL(("[%d] GRD Load Received objId %d  with load %f and toPE %d total_migrates %d total_migratesActual %d migrates_expected %d migrates_completed %d\n", CkMyPe(), objId, myStats->objData[objId].wallTime, fromPE, total_migrates, total_migratesActual, migrates_expected, migrates_completed));
}

void DiffusionLB::MigrationEnded() {
    // TODO: not deleted
    entered = true;
    DEBUGL(("[%d] GRD Migration Ended total_migrates %d total_migratesActual %d \n", CkMyPe(), total_migrates, total_migratesActual));
    msg = new(total_migrates,CkNumPes(),CkNumPes(),0) LBMigrateMsg;
    msg->n_moves = total_migrates;
    for(int i=0; i < total_migrates; i++) {
      MigrateInfo* item = (MigrateInfo*) migrateInfo[i];
      msg->moves[i] = *item;
      delete item;
      migrateInfo[i] = 0;
    }
    migrateInfo.clear();
    
    // Migrate messages from me to elsewhere
    for(int i=0; i < msg->n_moves; i++) {
        MigrateInfo& move = msg->moves[i];
        const int me = CkMyPe();
        if (move.from_pe == me && move.to_pe != me) {
	        lbmgr->Migrate(move.obj,move.to_pe);
        } else if (move.from_pe != me) {
	        DEBUGL(("[%d] error, strategy wants to move from %d to  %d\n",
		    me,move.from_pe,move.to_pe));
        }
    }
    if (CkMyPe() == 0) {
        CkCallback cb(CkIndex_DiffusionLB::MigrationDone(), thisProxy);
        CkStartQD(cb);
    }
}

//What does Cascading migrations do?
void DiffusionLB::CascadingMigration(LDObjHandle h, double load) {
    double threshold = THRESHOLD*avgLoadNeighbor/100.0;
    int minNode = -1;
    int myPos = neighborPos[CkNodeOfDiff(nodeFirst)];
    if(actualSend > 0) {
        double minLoad;
        // Send to max underloaded node
        for(int i = 0; i < neighborCount; i++) {
            if(balanced[i] == true && load <= toSendLoad[i] && (minNode == -1 || minLoad < toSendLoad[i])) {
                minNode = i;
                minLoad = toSendLoad[i];
            }
        }
        DEBUGL(("[%d] GRD Cascading Migration actualSend %d to node %d\n", CkMyPe(), actualSend, nbors[minNode]));
        if(minNode != -1 && minNode != myPos) {
            // Send load info to receiving load
            toSendLoad[minNode] -= load;
            if(toSendLoad[minNode] < threshold && balanced[minNode] == true) {
                balanced[minNode] = false;
                actualSend--; 
            }
            thisProxy[CkNodeFirstDiff(nbors[minNode])].LoadMetaInfo(h, load);
	        lbmgr->Migrate(h,CkNodeFirstDiff(nbors[minNode]));
        }
            
    }
    if(actualSend <= 0 || minNode == myPos || minNode == -1) {
        int minRank = -1;
        double minLoad = 0;
        for(int i = 0; i < nodeSize; i++) {
            if(minRank == -1 || loadPE[i] < minLoad) {
                minRank = i;
                minLoad = loadPE[i];
            }
        }
        DEBUGR(("[%d] GRD Cascading Migration actualSend %d sending to rank %d \n", CkMyPe(), actualSend, minRank));
        loadPE[minRank] += load;
        if(minRank > 0) {
	        lbmgr->Migrate(h, nodeFirst+minRank);
        }
    }
}

//What does this method do? - find out
void DiffusionLB::LoadMetaInfo(LDObjHandle h, double load) {
  //migrates_expected++;
  int idx = FindObjectHandle(h);
  if(idx == -1) {
    objectHandles.push_back(h);
    objectLoads.push_back(load);
  }
  else {
    CascadingMigration(h, load);
    objectHandles[idx] = objectHandles[objectHandles.size()-1];
    objectLoads[idx] = objectLoads[objectLoads.size()-1];
    objectHandles.pop_back();
    objectLoads.pop_back();
  }
}

void DiffusionLB::Migrated(LDObjHandle h, int waitBarrier)
{
  //migrates_completed++;
  if(CkMyPe() == nodeFirst) {
    thisProxy[CkMyPe()].MigratedHelper(h, waitBarrier);
  }
}

void DiffusionLB::MigratedHelper(LDObjHandle h, int waitBarrier) {
  DEBUGL(("[%d] GRD Migrated migrates_completed %d migrates_expected %d \n", CkMyPe(), migrates_completed, migrates_expected));
  int idx = FindObjectHandle(h);
  if(idx == -1) {
    objectHandles.push_back(h);
    objectLoads.push_back(-1);
  }
  else {
    CascadingMigration(h, objectLoads[idx]);
    objectHandles[idx] = objectHandles[objectHandles.size()-1];
    objectLoads[idx] = objectLoads[objectLoads.size()-1];
    objectHandles.pop_back();
    objectLoads.pop_back();
  }
}

void DiffusionLB::MigrationDone() {
    DEBUGR(("[%d] GRD Migration Done \n", CkMyPe()));
#if CMK_LBDB_ON
  migrates_completed = 0;
  total_migrates = 0;
  migrates_expected = -1;
  total_migratesActual = -1;
  avgLoadNeighbor = 0;
  if(CkMyPe() == 0) {
    end_lb_time = CkWallTimer();
    DEBUGL(("Strategy Time %f \n", end_lb_time - start_lb_time));
  }
  
  if(CkMyPe() == nodeFirst) {
    nodeStats->objData.clear();
    nodeStats->commData.clear();
    for(int i = 0; i < nodeSize; i++) {
      loadPE[i] = 0;
      numObjects[i] = 0;
      migratedTo[i] = 0;
      migratedFrom[i] = 0;
    }
    CallResumeClients();
  }

  // Increment to next step
  lbmgr->incStep();
  if(finalBalancing)
    lbmgr->ClearLoads();

  // if sync resume invoke a barrier
    if(!_lb_args.debug() || CkMyPe() != nodeFirst) {
    if (finalBalancing && _lb_args.syncResume()) {
    CkCallback cb(CkIndex_DiffusionLB::ResumeClients((CkReductionMsg*)(NULL)), 
        thisProxy);
    contribute(0, NULL, CkReduction::sum_int, cb);
    }
    else {
      thisProxy [CkMyPe()].ResumeClients(finalBalancing);
    }
  }
#endif
}

void DiffusionLB::ResumeClients(CkReductionMsg *msg) {
  ResumeClients(1);
  delete msg;
}

void DiffusionLB::CallResumeClients() {
  CmiAssert(_lb_args.debug());
//  DEBUGJ(("\n[PE-%d,Node-%d] GRD: Call Resume clients \n", CkMyPe(),CkMyNodeDiff()));
  CkCallback cbm(CkReductionTarget(DiffusionLB, SumIntBytes), thisProxy[0]);
//  contribute(sizeof(double), &internalBytes, CkReduction::sum_long, cbm);
  CkCallback cba(CkReductionTarget(DiffusionLB, SumExtBytes), thisProxy);
//  contribute(sizeof(double), &externalBytes, CkReduction::sum_long, cba);
  thisProxy[CkMyPe()].ResumeClients(finalBalancing);
}

void DiffusionLB::ResumeClients(int balancing) {
#if CMK_LBDB_ON

  if (CkMyPe() == 0 && balancing) {
    double end_lb_time = CkWallTimer();
    if (_lb_args.debug())
      DEBUGL(("%s> step %d finished at %f duration %f memory usage: %f\n",
          lbName(), step() - 1, end_lb_time, end_lb_time /*- strat_start_time*/,
          CmiMemoryUsage() / (1024.0 * 1024.0)));
  }

  lbmgr->ResumeClients();
#endif
}

// Aggregates the stats messages of PE into LDStats, Computes total load of node
void DiffusionLB::BuildStats()
{
#if DEBUG_K
    DEBUGL(("[%d] GRD Build Stats  and objects %lu\n", CkMyPe(), nodeStats->objData.size()));
#endif

    int n_objs = nodeStats->objData.size();
    int n_comm = nodeStats->commData.size();
//    DEBUGL(("\nn_objs=%d, n_comm=%d\n", n_objs,n_comm);
//    nodeStats->nprocs() = statsReceived;
    // allocate space
    nodeStats->objData.clear();
    nodeStats->from_proc.clear();
    nodeStats->to_proc.clear();
    nodeStats->commData.clear();
    int prev = 0;
    for(int i = 0; i < nodeSize; i++) {
      prefixObjects[i] = prev + numObjects[i];
      prev = prefixObjects[i];
    }

    nodeStats->objData.resize(n_objs);
    nodeStats->from_proc.resize(n_objs);
    nodeStats->to_proc.resize(n_objs);
    nodeStats->commData.resize(n_comm);
    objs.clear();
    objs.resize(n_objs);
       
    /*if(nodeKeys != NULL)
        delete[] nodeKeys;
    nodeKeys = new LDObjKey[nodeStats->n_objs];*/

    int nobj = 0;
    int ncom = 0;
    int nmigobj = 0;
    int start = nodeFirst;
    my_load = 0;
    my_loadAfterTransfer = 0;
    
    // copy all data in individual message to this big structure
    for (int pe=0; pe<statsReceived; pe++) {
      int i;
      CLBStatsMsg *msg = statsList[pe];
      if(msg == NULL) continue;
      for (i = 0; i < msg->objData.size(); i++) {
        nodeStats->from_proc[nobj] = nodeStats->to_proc[nobj] = start + pe;
        nodeStats->objData[nobj] = msg->objData[i];
        LDObjData &oData = nodeStats->objData[nobj];
#if MDEBUG
        DEBUGL(("\n[PE-%d]Adding vertex id %d", CkMyPe(), nobj));
#endif
        objs[nobj] = CkVertex(nobj, oData.wallTime, nodeStats->objData[nobj].migratable, nodeStats->from_proc[nobj]);
        my_load += msg->objData[i].wallTime;
        loadPE[pe] += msg->objData[i].wallTime;
        /*TODO Keys LDObjKey key;
        key.omID() = msg->objData[i].handle.omID;
        key.objID() =  msg->objData[i].handle.objID;
        nodeKeys[nobj] = key;*/
        if (msg->objData[i].migratable) 
            nmigobj++;
        nobj++;
      }
//        DEBUGL(("[%d] GRD BuildStats load of rank %d is %f \n", CkMyPe(), pe, loadPE[pe]);
      for (i = 0; i < msg->commData.size(); i++) {
        nodeStats->commData[ncom] = msg->commData[i];
        //nodeStats->commData[ncom].receiver.dest.destObj.destObjProc = msg->commData[i].receiver.dest.destObj.destObjProc; 
        int dest_pe = nodeStats->commData[ncom].receiver.lastKnown();
        //DEBUGL(("\n here dest_pe = %d\n", dest_pe);
        ncom++;
      }
      // free the memory TODO: Free the memory in Destructor
      delete msg;
      statsList[pe]=0;
    }
    my_loadAfterTransfer = my_load;
    nodeStats->n_migrateobjs = nmigobj;

    // Generate a hash with key object id, value index in objs vector
    nodeStats->deleteCommHash();
    nodeStats->makeCommHash();
}

void DiffusionLB::AddToList(CLBStatsMsg* m, int rank) {
  DEBUGR(("[%d] GRD Add To List num objects %d from rank %d load %f\n", CkMyPe(), m->n_objs, rank, m->total_walltime));
//    nodeStats->n_objs += m->n_objs;
//    nodeStats->n_comm += m->n_comm;
  nodeStats->objData.resize(nodeStats->objData.size()+m->objData.size());
  nodeStats->commData.resize(nodeStats->commData.size()+m->commData.size());
  numObjects[rank] = m->objData.size();
  statsList[rank] = m;
  
  struct ProcStats &procStat = nodeStats->procs[rank];
  procStat.pe = CkMyPe() + rank;	// real PE
  procStat.total_walltime = m->total_walltime;
  procStat.idletime = m->idletime;
  procStat.bg_walltime = m->bg_walltime;
  #if CMK_LB_CPUTIMER
  procStat.total_cputime = m->total_cputime;
  procStat.bg_cputime = m->bg_cputime;
  #endif
  procStat.pe_speed = m->pe_speed;		// important
  procStat.available = true;
  procStat.n_objs = m->objData.size();
}
#include "DiffusionLB.def.h"

