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

#include "ck.h"
#include "ckgraph.h"
#include "envelope.h"
//#include "LBDBManager.h"
//#include "LBSimulation.h"
#include "elements.h"
#include "Heap_helper.C"
#define DEBUGF(x) CmiPrintf x;
#define DEBUGR(x) /*CmiPrintf x*/;
#define DEBUGL(x) /*CmiPrintf x*/;
#define NUM_NEIGHBORS 2
#define ITERATIONS 4

#define NODESIZE 2

// Percentage of error acceptable.
#define THRESHOLD 2

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

using std::vector;
// TODO'S: 
// Topology
// Non migratable objects

void DiffusionLB::staticMigrated(void* data, LDObjHandle h, int waitBarrier)
{
  DiffusionLB *me = (DiffusionLB*)(data);
  me->Migrated(h, waitBarrier);
}

void DiffusionLB::staticAtSync(void* data)
{
  DiffusionLB *me = (DiffusionLB*)(data);
//  me->AtSync();
}

void DiffusionLB::Strategy(const DistBaseLB::LDStats* const stats) {
  if (CkMyPe() == 0 && _lb_args.debug() >= 1) {
    double start_time = CmiWallTimer();
    CkPrintf("In DiffusionLB strategy at %lf\n", start_time);
  }
  lb_started = false;
  AtSync();
}
// preprocess topology information: only done once
void DiffusionLB::preprocess(const int explore_limit)
{
}


DiffusionLB::DiffusionLB(CkMigrateMessage *m) : CBase_DiffusionLB(m) {
}

DiffusionLB::DiffusionLB(const CkLBOptions &opt) : CBase_DiffusionLB(opt) {
#if CMK_LBDB_ON
  lbname = "DiffusionLB";
  if (CkMyPe() == 0)
      CkPrintf("[%d] Diffusion created\n",CkMyPe());
  if (_lb_args.statsOn()) lbmgr->CollectStatsOn();
  InitLB(opt);
#endif
}

DiffusionLB::~DiffusionLB()
{
#if CMK_LBDB_ON
  delete [] statsList;
  delete nodeStats;
  //delete [] myStats->objData;
  //delete [] myStats->commData;
  delete myStats;
  delete[] gain_val;
  delete[] obj_arr;
  lbmgr = CProxy_LBManager(_lbmgr).ckLocalBranch();
  if (lbmgr)
    lbmgr->RemoveStartLBFn(startLbFnHdl); 
#endif
}

void DiffusionLB::InitLB(const CkLBOptions &opt) {
#if DEBUG_K
  DEBUGF(("[%d] InitLB\n", CkMyPe()));
#endif
  thisProxy = CProxy_DiffusionLB(thisgroup);
  numNodes = CkNumNodesDiff();
//    receiver = lbmgr->AddLocalBarrierReceiver((LDBarrierFn)(staticAtSync),
//                (void*)(this));
//    notifier = lbmgr->NotifyMigrated((LDMigratedFn)(staticMigrated), (void*)(this));
  startLbFnHdl = lbmgr->AddStartLBFn(this, &DiffusionLB::ProcessAtSync);

  myspeed = lbmgr->ProcessorSpeed();
  // TODO: Initialize all class variables
  loadReceived = 0;
  statsReceived = 0;
  total_migrates = 0;
  total_migratesActual = -1;
  migrates_expected = -1;
  migrates_completed = 0;
  myStats = new DistBaseLB::LDStats;
  //myStats->objData = NULL;
  //myStats->commData = NULL;
  nodeFirst = CkNodeFirstDiff(CkMyNodeDiff());
  notif = 0;
#if DEBUG_K
  CkPrintf("\n[PE-%d] nodeFirst = %d", CkMyPe(), nodeFirst);
#endif
  ComputeNeighbors();
  if(CkMyPe() == nodeFirst) {
    gain_val = NULL;
    obj_arr = NULL;
    my_load = 0;
    my_loadAfterTransfer = 0;
    toSend = 0;
    statsList = new CLBStatsMsg*[nodeSize];
    nodeStats = new BaseLB::LDStats(nodeSize);
    loadPE.reserve(nodeSize);
    loadPEBefore.reserve(nodeSize);
    numObjects.reserve(nodeSize);
    loadNeighbors.reserve(neighborCount);
    prefixObjects.reserve(nodeSize);
    migratedTo.reserve(nodeSize);
    migratedFrom.reserve(nodeSize);
    for(int i = 0; i < nodeSize; i++) {
      loadPE[i] = 0;
      loadPEBefore[i] = 0;
      numObjects[i] = 0;
      migratedTo[i] = 0;
      migratedFrom[i] = 0;
      prefixObjects[i] = 0;
    }
  }
}

void DiffusionLB::AtSync() {
#if DEBUG_K
  CkPrintf("\n[PE-%d]In DiffusionLB::AtSync()\n", CkMyPe()); fflush(stdout);
#endif
#if CMK_LBDB_ON
  if (!QueryBalanceNow(step()) || CkNumPes() == 1) {
    finalBalancing = 0;
    MigrationDone();
    return;
  }
  finalBalancing = 1;
  if(CkMyPe() == 0 && _lb_args.debug()) {
    maxB = 0;
    minB = -1;
    avgB = 0;
    maxA = 0;
    minA = -1;
    avgA = 0;
    receivedStats = 0;
    internalBeforeFinal = 0;
    externalBeforeFinal = 0;
    internalAfterFinal = 0;
    externalAfterFinal = 0;
  }
  internalBefore = 0;
  externalBefore = 0;
  internalAfter = 0;
  externalAfter = 0;
  migrates = 0;
  migratesNode = 0;
  objectHandles.clear();
  objectLoads.clear();
  // TODO: Check is it is the first load balancing step and then only 
  // perform this sending and QD
  if(step() == 0) {
    sendToNeighbors.clear();
    //toSend = 0;
//    CkPrintf("\nCkMyPe()=%d, nodeFirst=%d", CkMyPe(), nodeFirst);
    CkCallback cb(CkIndex_DiffusionLB::PEStarted(), thisProxy[0]);
    contribute(cb); 
  }
  else {
    //if(CmiNodeAlive(CkMyPe())){
        thisProxy[CkMyPe()].ProcessAtSync();
    //}
  }
#endif
}

void DiffusionLB::PEStarted() {
  if(CkMyPe() == 0) {
    CkCallback cb(CkIndex_DiffusionLB::ProcessAtSync(), thisProxy);
    CkStartQD(cb); 
  }
}

void DiffusionLB::AddNeighbor(int node) {
  toSend++;
#if 1//DEBUG_K
  CkPrintf("[PE-%d, Node-%d] My %luth neighbor node is %d\n", CkMyPe(), CkMyNodeDiff(), sendToNeighbors.size(), node);
#endif
  if(sendToNeighbors.size() > neighborCount)
    CkPrintf("\n[PE-%d,node-%d]Adding nbors (node-%d) beyond count!! %lu>(of expected count %d)\n", CkMyPe(), CkMyNodeDiff(), node, sendToNeighbors.size(), neighborCount);
  sendToNeighbors.push_back(node);
}

void DiffusionLB::ProcessAtSync()
{
#if CMK_LBDB_ON
  start_lb_time = 0;
#if DEBUG_K
  CkPrintf("[%d] ProcessAtSync()", CkMyPe());
#endif

  if(step() == 0) { 
//    toSendLoad.resize(neighborCount);
//    toReceiveLoad.resize(sendToNeighbors.size());

//    for(int i = 0; i < sendToNeighbors.size(); i++) {
      //neighborPosReceive[sendToNeighbors[i]] = i;
//      toReceiveLoad[i] = 0;
//    }
  }
//  for(int i = 0; i < neighborCount; i++)
//    toSendLoad[i] = 0;
//  for(int i = 0; i < sendToNeighbors.size(); i++)
//    toReceiveLoad[i] = 0;

  if (CkMyPe() == 0) {
    start_lb_time = CkWallTimer();
    if (_lb_args.debug())
      CkPrintf("[%s] Load balancing step -%d starting at %f\n",
                              lbName(), step(), CkWallTimer());
  }

#if DEBUG_K
  CkPrintf("\nBefore AssembleStats()");
#endif

  // assemble LB database
  statsmsg = AssembleStats();
  if(statsmsg == NULL)
    CkPrintf("!!!Its null!!!\n");

 // CkMarshalledCLBStatsMessage *marshmsg;
//  marshmsg = new CkMarshalledCLBStatsMessage(statsmsg);
 
    // send to parent
  CkPrintf("[%d] Sending to parent #%d ReceiveStats\n", CkMyPe(), nodeFirst); fflush(stdout);
  CkMarshalledCLBStatsMessage marshmsg(statsmsg);
//  marshmsg = new CkMarshalledCLBStatsMessage(statsmsg);
  thisProxy[nodeFirst].ReceiveStats(marshmsg);
  if(CkMyPe() != nodeFirst) {
    CkCallback cb(CkIndex_DiffusionLB::createNeighbors(), thisProxy);
    contribute(cb);
  }
}

void DiffusionLB::doneNborExng() {
  if(CkMyPe() == CkNodeFirstDiff(CkMyNodeDiff()) && step() == 0) {
      loadNeighbors.clear();
      toReceiveLoad.clear();
      toSendLoad.clear();
      sendToNeighbors.clear();
      loadNeighbors.reserve(neighborCount);
      toReceiveLoad.resize(neighborCount);
      toSendLoad.resize(neighborCount);
      sendToNeighbors.reserve(neighborCount);
      CkPrintf("\n[PE-%d] setting sendToNeighbors size to %d", CkMyPe(), neighborCount);
      std::string nbor_nodes = " ";
      for(int i = 0; i < neighborCount; i++) {
        //Add your neighbors node-id as your neighbor
         nbor_nodes += "node-"+ std::to_string(nbors[i])+", ";
        AddNeighbor(nbors[i]);
      }
      CkPrintf("\n[PE-%d,Node-%d] my neighbors: %s (#%d neighbors)\n", CkMyPe(), CkMyNodeDiff(), nbor_nodes.c_str(), neighborCount);
  }

  // send to parent
  CkPrintf("[%d] Sending to parent #%d ReceiveStats\n", CkMyPe(), nodeFirst); fflush(stdout);
//  CkMarshalledCLBStatsMessage marshmsg(statsmsg);
//  marshmsg = new CkMarshalledCLBStatsMessage(statsmsg);
//  thisProxy[nodeFirst].ReceiveStats(*marshmsg);
  if(CkMyPe() == CkNodeFirstDiff(CkMyNodeDiff()))
  thisProxy[CkMyPe()].iterate();
#endif
}

void DiffusionLB::ComputeNeighbors() {
  //TODO: Use application graph for computing node neighbors
  //Assuming nodes are neighbors in a line
  nodeSize = CkMyNodeSizeDiff(); 
  nodeFirst = CkMyNodeDiff()*nodeSize;
  // TODO: Juan's topology aware mapping
  neighborCount = NUM_NEIGHBORS/2;
/*
  neighbors.resize(neighborCount);
  neighbors[0] = CkMyNodeDiff()-1;
  if(neighbors[0] < 0) neighbors[0] = CkNumNodesDiff()-1;
  neighbors[1] = (CkMyNodeDiff()+1)%CkNumNodesDiff();
#if DEBUG_K
  for(int i=0;i<neighbors.size();i++)
    CkPrintf("\nPE-%d, neighbor node = %d", CkMyPe(), neighbors[i]);
#endif
*/    
#if 0
    preprocess(NUM_NEIGHBORS);

    int pdims[4];
    TopoManager_getPeCoordinates(CkMyPe(), pdims);
    std::vector<TCoord*> closest = closest_coords[TCoord(pdims).idx]; 
    TCoord &c = coord_table[TCoord(pdims).idx];
    nodeSize = c.ppn;
    nodeFirst = c.p;
    DEBUGR(("[%d] nodeSize %d and nodeFirst %d \n", CkMyPe(), nodeSize, nodeFirst)); 
    int dist = 0;
    if(CkMyPe() == nodeFirst)
    for(int i = 0; i < closest.size(); i++) {
        if(peNodes.find(closest[i]->p) == peNodes.end())
            CkAbort("all rank 0 pe's not included\n");
        int node = peNodes[closest[i]->p];
        if(neighborPos.find(node) == neighborPos.end()) {
            neighbors.push_back(node);
            neighborPos[neighbors[dist]] = dist;
            dist++;
            CkPrintf("[%d] GRD: neighbor is %d \n", CkMyPe(), node); 
        }
    }
    neighborCount = dist;
#endif
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
    CkPrintf("\nPE-%d Error!!!!!", CkMyPe());
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

//    myStats->n_objs = lbmgr->GetObjDataSz();
//    if (myStats->objData != NULL) { 
#if DEBUG_K
  CkPrintf("Freeing \n");
#endif
//    myStats->objData.erase();//delete[] myStats->objData;
//    myStats->objData = NULL;}
  myStats->objData.resize(lbmgr->GetObjDataSz());// = new LDObjData[myStats->n_objs];
  lbmgr->GetObjData(myStats->objData.data());

//    myStats->n_comm = lbmgr->GetCommDataSz();
//    delete[] myStats->commData;
  myStats->commData.resize(lbmgr->GetCommDataSz());// = new LDCommData[myStats->n_comm];
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
//  msg->pe_speed = myspeed;
  // number of pes
  statsMsg->pe_speed = myStats->pe_speed;

  //statsMsg->n_objs = osz;
  lbmgr->GetObjData(statsMsg->objData.data());
  //statsMsg->n_comm = csz;
  lbmgr->GetCommData(statsMsg->commData.data());

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
  DEBUGF(("[%d] GRD ReceiveStats from pe %d\n", CkMyPe(), m->from_pe));
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
    CkCallback cb(CkIndex_DiffusionLB::createNeighbors(), thisProxy);
    contribute(cb);
  }
}
void DiffusionLB::createNeighbors(){
    if(step() == 0 && CkMyPe()==CkNodeFirstDiff(CkMyNodeDiff())) {
    long ebytes[CkNumNodesDiff()];
    std::fill_n(ebytes, CkNumNodesDiff(), 0);
    nbors = new int[NUM_NEIGHBORS+CkNumNodesDiff()];
    for(int i=0;i<CkNumNodesDiff();i++)
      nbors[i] = -1;
    neighborCount = NUM_NEIGHBORS/2;
    CkPrintf("\nedges = %lu", nodeStats->commData.size());
    for(int edge = 0; edge < nodeStats->commData.size(); edge++) {
      LDCommData &commData = nodeStats->commData[edge];
      if( (!commData.from_proc()) && (commData.recv_type()==LD_OBJ_MSG) ) {
        LDObjKey from = commData.sender;
        LDObjKey to = commData.receiver.get_destObj();

        // Check the possible values of lastKnown.
        int toPE = commData.receiver.lastKnown();
        int toNode = CkNodeOfDiff(toPE);
        CkPrintf("\ntoPE = %d",toPE);
        if(CkMyNodeDiff() != toNode && toNode!= -1) {
          ebytes[toNode] += commData.bytes;
        }
      }
    }
    for(int i=0;i<CkNumNodesDiff();i++)
      CkPrintf("\n[PE-%d,Node-%d] ebytes[to node %d] = %lu", CkMyPe(), CkMyNodeDiff(), i, ebytes[i]);
    sortArr(ebytes, CkNumNodesDiff(), nbors);
    CkPrintf("\n[PE-%d, node-%d], my largest comm neighbors are %d,%d\n", CkMyPe(), CkMyNodeDiff(), nbors[0], nbors[1]);
  }

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
        CkPrintf("\n[PE-%d], notifying node %d [PE-%d]\n", CkMyPe(), i, CkNodeFirstDiff(i));
  #endif
      thisProxy[CkNodeFirstDiff(i)].notifyNeighbor(isNbor, CkMyNodeDiff());
    }
  } else
    doneNborExng();

    statsReceived = 0;
//    thisProxy[CkMyPe()].iterate();

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
/*
void Diffusion::ReceiveLoadInfo(double load, int node) {
    DEBUGR(("[%d] GRD Receive load info, load %f node %d loadReceived %d neighborCount %d\n", CkMyPe(), load, node, loadReceived, neighborCount));
    int pos = neighborPos[node];
    loadNeighbors[pos] = load;
    loadReceived++;

    if(loadReceived == neighborCount) {
        loadReceived = 0;     
        avgLoadNeighbor = average();
        DEBUGR(("[%d] GRD Received all loads of node, avg is %f and my_load %f \n", CkMyPe(), avgLoadNeighbor, my_load));
        double threshold = THRESHOLD*avgLoadNeighbor/100.0;
        if(my_load > avgLoadNeighbor + threshold) {
            LoadBalancing();
        }
        else if (CkMyPe() == 0) {
            CkCallback cb(CkIndex_Diffusion::DoneNodeLB(), thisProxy);
            CkStartQD(cb);
        }
    }
}
*/
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
    int node = nbors[i];
    if(neighborPosReceive.find(node) != neighborPosReceive.end()) {
      // One of them will become negative
      int pos = neighborPosReceive[node];
      toSendLoad[i] -= toReceiveLoad[pos];
      if(toSendLoad[i] > 0)
          res= true;
      toReceiveLoad[pos] -= toSendLoad[i];
    }
    CkPrintf("[%d] Diff: To Send load to node %d load %f res %d\n", CkMyPe(), node, toSendLoad[i], res);
    CkPrintf("[%d] Diff: To Send load to node %d load %f res %d myLoadB %f\n", CkMyPe(), node, toSendLoad[i], res, my_loadAfterTransfer);
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

void DiffusionLB::PseudoLoadBalancing() {
  CkPrintf("[PE-%d] Pseudo Load Balancing , iteration %d my_load %f my_loadAfterTransfer %f avgLoadNeighbor %f\n", CkMyPe(), itr, my_load, my_loadAfterTransfer, avgLoadNeighbor);
  double threshold = THRESHOLD*avgLoadNeighbor/100.0;
  
  double totalOverload = my_load - avgLoadNeighbor;
  avgLoadNeighbor = (avgLoadNeighbor+my_load)/2;
  double totalUnderLoad = 0.0;
  double thisIterToSend[neighborCount];
  for(int i = 0 ;i < neighborCount; i++)
    thisIterToSend[i] = 0.;
  if(totalOverload > 0)
    for(int i = 0; i < neighborCount; i++) {
      thisIterToSend[i] = 0;
      if(loadNeighbors[i] < (avgLoadNeighbor - threshold)) {
        thisIterToSend[i] = avgLoadNeighbor - loadNeighbors[i];
        totalUnderLoad += avgLoadNeighbor - loadNeighbors[i];
        CkPrintf("[PE-%d] iteration %d thisIterToSend %f avgLoadNeighbor %f loadNeighbors[i] %f to node %d\n",
                CkMyPe(), itr, thisIterToSend[i], avgLoadNeighbor, loadNeighbors[i], nbors[i]);
      }
    }
  if(totalUnderLoad > 0 && totalOverload > 0 && totalUnderLoad > totalOverload)
    totalOverload += threshold;
  else
    totalOverload = totalUnderLoad;
  DEBUGL(("[%d] GRD: Pseudo Load Balancing Sending, iteration %d totalUndeload %f totalOverLoad %f my_loadAfterTransfer %f\n", CkMyPe(), itr, totalUnderLoad, totalOverload, my_loadAfterTransfer));
  for(int i = 0; i < neighborCount; i++) {
    if(totalOverload > 0 && totalUnderLoad > 0 && thisIterToSend[i] > 0) {
      CkPrintf("[%d] GRD: Pseudo Load Balancing Sending, iteration %d node %d(pe-%d) toSend %lf totalToSend %lf\n", CkMyPe(), itr, nbors[i], CkNodeFirstDiff(nbors[i]), thisIterToSend[i], (thisIterToSend[i]*totalOverload)/totalUnderLoad);
      thisIterToSend[i] *= totalOverload/totalUnderLoad;
      toSendLoad[i] += thisIterToSend[i];
    }
    if(my_load - thisIterToSend[i] < 0)
      CkAbort("Get out");
    my_load -= thisIterToSend[i];
    thisProxy[CkNodeFirstDiff(nbors[i])].PseudoLoad(itr, thisIterToSend[i], CkMyNodeDiff());
  }
}

int DiffusionLB::findNborIdx(int node) {
  CkPrintf("\n[PE-%d]Looking for node %d in nbor array of size %lu", CkMyPe(), node, sendToNeighbors.size());
  for(int i=0;i<sendToNeighbors.size();i++)
    if(sendToNeighbors[i] == node)
      return i;
  return -1;
}

#define SELF_IDX NUM_NEIGHBORS
#define EXT_IDX NUM_NEIGHBORS+1
void DiffusionLB::LoadBalancing() {
  int n_objs = nodeStats->objData.size();
  CkPrintf("[%d] GRD: Load Balancing w objects size = %d \n", CkMyPe(), n_objs);

//  Iterate over the comm data and for each object, store its comm bytes
//  to other neighbor nodes and own node.

  //objectComms maintains the comm bytes for each object on this node
  //with the neighboring node
  //we also maintain comm within this node and comm bytes outside
  //(of this node and neighboring nodes)
  vector<vector<int>> objectComms(n_objs);
//  objectComms.resize(n_objs);

  if(gain_val != NULL)
      delete[] gain_val;
  gain_val = new int[n_objs];
  memset(gain_val, -1, n_objs);

  CkPrintf("\n[PE-%d] n_objs=%d", CkMyPe(), n_objs);

  for(int i = 0; i < n_objs; i++) {
//    CkPrintf("\nPE-%d objid= %" PRIu64 ", vrtx id=%d", CkMyPe(), nodeStats->objData[i].objID(), objs[i].getVertexId());
    objectComms[i].resize(neighborCount+1);
    for(int j = 0; j < neighborCount+1; j++)
      objectComms[i][j] = 0;
  }

  // TODO: Set objectComms to zero initially
  int obj = 0;
  for(int edge = 0; edge < nodeStats->commData.size(); edge++) {
    LDCommData &commData = nodeStats->commData[edge];
    // ensure that the message is not from a processor but from an object
    // and that the type is an object to object message
    if( (!commData.from_proc()) && (commData.recv_type()==LD_OBJ_MSG) ) {
      LDObjKey from = commData.sender;
      LDObjKey to = commData.receiver.get_destObj();
      int fromNode = CkMyNodeDiff();//peNodes[nodeFirst]; //Originating from my node? - q

      // Check the possible values of lastKnown.
      int toPE = commData.receiver.lastKnown();
      int toNode = CkNodeOfDiff(toPE);
      //store internal bytes in the last index pos ? -q
      if(fromNode == toNode) {
//        int pos = neighborPos[toNode];
        int nborIdx = SELF_IDX;// why self id?
        int fromObj = nodeStats->getHash(from);
        int toObj = nodeStats->getHash(to);
        //DEBUGR(("[%d] GRD Load Balancing from obj %d and to obj %d and total objects %d\n", CkMyPe(), fromObj, toObj, nodeStats->n_objs));
        objectComms[fromObj][nborIdx] += commData.bytes;
        // lastKnown PE value can be wrong.
        if(toObj != -1) {
          objectComms[toObj][nborIdx] += commData.bytes; 
          internalBefore += commData.bytes;
        }
        else
          externalBefore += commData.bytes;
      }
      else { // External communication? - q
        externalBefore += commData.bytes;
        int nborIdx = findNborIdx(toNode);
        if(nborIdx == -1)
          nborIdx = EXT_IDX;//Store in last index if it is external bytes going to
//        non-immediate neighbors? -q
        if(fromNode == CkMyNodeDiff()/*peNodes[nodeFirst]*/) {//ensure bytes are going from my node? -q
          int fromObj = nodeStats->getHash(from);
          CkPrintf("[%d] GRD Load Balancing from obj %d and pos %d\n", CkMyPe(), fromObj, nborIdx);
          objectComms[fromObj][nborIdx] += commData.bytes;
          obj++;
        }
      }
    }
    } // end for

    // calculate the gain value, initialize the heap.
    internalAfter = internalBefore;
    externalAfter = externalBefore;
    double threshold = THRESHOLD*avgLoadNeighbor/100.0;

    actualSend = 0;
    balanced.resize(toSendLoad.size());
    CkPrintf("\nIterating through toSendLoad of size %lu", toSendLoad.size());
    for(int i = 0; i < toSendLoad.size(); i++) {
      balanced[i] = false;
      if(toSendLoad[i] > 0) {
        balanced[i] = true;
        actualSend++;
      }
    }
    CkPrintf("\n[PE-%d,Node-%d]actualSend = %d", CkMyPe(), CkMyNodeDiff(), actualSend);

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
      CkPrintf("\n[PE-%d] my_loadAfterTransfer = %lf, actualSend=%d\n", CkMyPe(),my_loadAfterTransfer,actualSend);

      //return;
      while(my_loadAfterTransfer > 0 && actualSend > 0) {
        counter++;
        //pop the object id with the least gain (i.e least internal comm compared to ext comm)

        v_id = heap_pop(obj_heap, ObjCompareOperator(&objs, gain_val), heap_pos);
        
        CkPrintf("\n On PE-%d, popped v_id = %d", CkMyPe(), v_id);
   
        /*If the heap becomes empty*/
        if(v_id==-1)          
            break;
        double currLoad = objs_cpy[v_id].getVertexLoad();
#if 0
        if(!objs[v_id].isMigratable()) {
          CkPrintf("not migratable \n");
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
        CkPrintf("\n[PE-%d] maxi = %d", CkMyPe(), maxi);
          
        if(maxi != -1) {
#if 1
          migrates++;
          int pos = neighborPos[CkNodeOfDiff(nodeFirst)];
          internalAfter -= comm[pos];
          internalAfter += comm[maxi];
          externalAfter += comm[pos];
          externalAfter -= comm[maxi];
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
              CkPrintf("\n%d!=%d", objId, v_id);fflush(stdout);
              CmiAbort("objectIds dont match \n");
          }
          int pe = GetPENumber(objId);
          migratedFrom[pe]++;
          int initPE = CkNodeFirstDiff(CkMyNodeDiff()) + pe;
          loadPE[pe] -= currLoad;
          numObjects[pe]--;
          CkPrintf("[%d] GRD: Load Balancing object load %f to node %d and from pe %d and objID %d\n", CkMyPe(), currLoad, node, initPE, objId);
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
          CkPrintf("[%d] maxi is negative currLoad %f \n", CkMyPe(), currLoad);
        } 
      } //end of while
      CkPrintf("[%d] GRD: Load Balancing total load sent during LoadBalancing %f actualSend %d myloadB %f v_id %d counter %d nobjs %lu \n",
          CkMyPe(), totalSent, actualSend, my_loadAfterTransfer, v_id, counter, nodeStats->objData.size());
      for (int i = 0; i < neighborCount; i++) {
        CkPrintf("[%d] GRD: Load Balancing total load sent during LoadBalancing toSendLoad %f node %d\n", CkMyPe(), toSendLoad[i], nbors[i]);
        }
      }//end of if
      // TODO: Put QD in intra node
      /* Start quiescence detection at PE 0.
      if (CkMyPe() == 0) {
          CkCallback cb(CkIndex_Diffusion::DoneNodeLB(), thisProxy);
          CkStartQD(cb);
      }*/
}

// Load is sent from overloaded to underloaded nodes, now we should load balance the PE's within the node
void DiffusionLB::DoneNodeLB() {
  entered = false;
  if(CkMyPe() == nodeFirst) {
    CkPrintf("[%d] GRD: DoneNodeLB \n", CkMyPe());
    double avgPE = averagePE();

    // Create a max heap and min heap for pe loads
    vector<double> objectSizes;
    vector<int> objectIds;
    minHeap minPes(nodeSize);
    double threshold = THRESHOLD*avgPE/100.0;
    
    for(int i = 0; i < nodeSize; i++) {
      if(loadPE[i] > avgPE + threshold) {
        CkPrintf("[%d] GRD: DoneNodeLB rank %d is overloaded with load %f\n", CkMyPe(), i, loadPE[i]);
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
        CkPrintf("[%d] GRD: DoneNodeLB rank %d is underloaded with load %f\n", CkMyPe(), i, loadPE[i]);
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
  CkPrintf("[%d] GRD Load Received objId %d  with load %f and toPE %d total_migrates %d total_migratesActual %d migrates_expected %d migrates_completed %d\n", CkMyPe(), objId, myStats->objData[objId].wallTime, fromPE, total_migrates, total_migratesActual, migrates_expected, migrates_completed);
}

void DiffusionLB::MigrationEnded() {
    // TODO: not deleted
    entered = true;
    CkPrintf("[%d] GRD Migration Ended total_migrates %d total_migratesActual %d \n", CkMyPe(), total_migrates, total_migratesActual);
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
	        CkPrintf("[%d] error, strategy wants to move from %d to  %d\n",
		    me,move.from_pe,move.to_pe);
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
        CkPrintf("[%d] GRD Cascading Migration actualSend %d to node %d\n", CkMyPe(), actualSend, nbors[minNode]);
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
    CkPrintf("[%d] GRD Migrated migrates_completed %d migrates_expected %d \n", CkMyPe(), migrates_completed, migrates_expected);
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

void DiffusionLB::PrintDebugMessage(int len, double* result) {
    avgB += result[2];
    if(result[0] > maxB) {
        maxB=result[0];
        maxPEB = (int)result[12];
    }
    if(minB == -1 || result[1] < minB) {
        minB = result[1];
        minPEB = (int)result[11];
    }
    avgA += result[5];
    if(result[3] > maxA) {
        maxA=result[3];
        maxPEA = (int)result[14];
    }
    if(minA == -1 || result[4] < minA) {
        minA = result[4];
        minPEA = (int)result[13];
    }
    internalBeforeFinal += result[6];
    externalBeforeFinal += result[7];
    internalAfterFinal += result[8];
    externalAfterFinal += result[9];
    migrates += result[10];
    
    receivedStats++;
    CkPrintf("\nPrints on PE %d, receivedStats = %d\n", CkMyPe(), receivedStats);
    if(receivedStats == numNodes) {
        receivedStats = 0;
        avgB = avgB /CkNumPes();
        avgA = avgA / CkNumPes();
        CkPrintf("Max PE load before %f(%d), after %f(%d) \n", maxB, maxPEB, maxA, maxPEA);
        CkPrintf("Min PE load before %f(%d), after %f(%d) \n", minB, minPEB, minA, minPEA);
        CkPrintf("Avg PE load before %f, after %f \n", avgB, avgA);
        CkPrintf("Internal Communication before %f, after %f \n", internalBeforeFinal, internalAfterFinal);
        CkPrintf("External communication before %f, after %f \n", externalBeforeFinal, externalAfterFinal);
        CkPrintf("Number of migrations across nodes %d \n", migrates);
        for(int i = 0; i < numNodes; i++) {
          CkPrintf("\nnodes[%d] = %d",i, i);
          thisProxy[CkNodeFirstDiff(i)].CallResumeClients();
        }
    }
    fflush(stdout);
}

void DiffusionLB::MigrationDone() {
    DEBUGR(("[%d] GRD Migration Done \n", CkMyPe()));
#if CMK_LBDB_ON
  migrates_completed = 0;
  total_migrates = 0;
  migrates_expected = -1;
  total_migratesActual = -1;
  avgLoadNeighbor = 0;
  //myStats->objData = NULL;
//  myStats->objData.erase();
//  delete[] myStats->objData;
//  myStats->commData.erase();
//  delete[] myStats->commData;
//  myStats->commData = NULL; 
    if(CkMyPe() == 0) {
        end_lb_time = CkWallTimer();
        CkPrintf("Strategy Time %f \n", end_lb_time - start_lb_time);
    }
    
    if(CkMyPe() == nodeFirst) {
        double minLoadB = loadPEBefore[0];
        double maxLoadB = loadPEBefore[0];
        double sumBefore = 0.0;
        double minLoadA = loadPE[0];
        double maxLoadA = loadPE[0];
        double sumAfter = 0.0;
        double maxPEB = nodeFirst;
        double maxPEA = nodeFirst;
        double minPEB = nodeFirst;
        double minPEA = nodeFirst;
        if (_lb_args.debug()) {
            for(int i = 0; i < nodeSize; i++) {
                CkPrintf("[%d] GRD: load of PE before: %f after: %f\n",CkMyPe()+i, loadPEBefore[i], loadPE[i] );
                if(minLoadB > loadPEBefore[i]) {
                    minLoadB = loadPEBefore[i];
                    minPEB = nodeFirst + i;
                }
                if(maxLoadB < loadPEBefore[i]) {
                    maxLoadB = loadPEBefore[i];
                    maxPEB = nodeFirst+i;
                }
                sumBefore += loadPEBefore[i];
                if(minLoadA > loadPE[i]) {
                    minLoadA = loadPE[i];
                    minPEA = nodeFirst + i;
                }
                if(maxLoadA < loadPE[i]) {
                    maxLoadA = loadPE[i];
                    maxPEA = nodeFirst + i;
                }
                sumAfter += loadPE[i];
            }
            double loads[15];
            loads[0] = maxLoadB;
            loads[1] = minLoadB;
            loads[2] = sumBefore;
            loads[3] = maxLoadA;
            loads[4] = minLoadA;
            loads[5] = sumAfter;
            loads[6] = internalBefore;
            loads[7] = externalBefore;
            loads[8] = internalAfter;
            loads[9] = externalAfter;
            loads[10] = migrates;
            loads[11] = minPEB;
            loads[12] = maxPEB;
            loads[13] = minPEA;
            loads[14] = maxPEA;
            DEBUGL(("[%d] PE's %f %f %f %f \n", CkMyPe(), minPEB, maxPEB, minPEA, maxPEA));
            thisProxy[0].PrintDebugMessage(15, loads);
        }
    
        nodeStats->objData.clear();
        nodeStats->commData.clear();
        for(int i = 0; i < nodeSize; i++) {
            loadPE[i] = 0;
            loadPEBefore[i] = 0;
            numObjects[i] = 0;
            migratedTo[i] = 0;
            migratedFrom[i] = 0;
        }
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
    CkPrintf("\n[PE-%d,Node-%d] GRD: Call Resume clients \n", CkMyPe(),CkMyNodeDiff());
    thisProxy[CkMyPe()].ResumeClients(finalBalancing);
}

void DiffusionLB::ResumeClients(int balancing) {
#if CMK_LBDB_ON

  if (CkMyPe() == 0 && balancing) {
    double end_lb_time = CkWallTimer();
    if (_lb_args.debug())
      CkPrintf("%s> step %d finished at %f duration %f memory usage: %f\n",
          lbName(), step() - 1, end_lb_time, end_lb_time /*- strat_start_time*/,
          CmiMemoryUsage() / (1024.0 * 1024.0));
  }

  lbmgr->ResumeClients();
#endif
}

// Aggregates the stats messages of PE into LDStats, Computes total load of node
void DiffusionLB::BuildStats()
{
#if DEBUG_K
    CkPrintf("[%d] GRD Build Stats  and objects %lu\n", CkMyPe(), nodeStats->objData.size());
#endif

    int n_objs = nodeStats->objData.size();
    int n_comm = nodeStats->commData.size();
//    CkPrintf("\nn_objs=%d, n_comm=%d\n", n_objs,n_comm);
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
        CkPrintf("\n[PE-%d]Adding vertex id %d", CkMyPe(), nobj);
#endif
        objs[nobj] = CkVertex(nobj, oData.wallTime, nodeStats->objData[nobj].migratable, nodeStats->from_proc[nobj]);
        my_load += msg->objData[i].wallTime;
        loadPE[pe] += msg->objData[i].wallTime;
        loadPEBefore[pe] += msg->objData[i].wallTime;
        /*TODO Keys LDObjKey key;
        key.omID() = msg->objData[i].handle.omID;
        key.objID() =  msg->objData[i].handle.objID;
        nodeKeys[nobj] = key;*/
        if (msg->objData[i].migratable) 
            nmigobj++;
        nobj++;
      }
//        CkPrintf("[%d] GRD BuildStats load of rank %d is %f \n", CkMyPe(), pe, loadPE[pe]);
        for (i = 0; i < msg->commData.size(); i++) {
            nodeStats->commData[ncom] = msg->commData[i];
            //nodeStats->commData[ncom].receiver.dest.destObj.destObjProc = msg->commData[i].receiver.dest.destObj.destObjProc; 
            int dest_pe = nodeStats->commData[ncom].receiver.lastKnown();
            //CkPrintf("\n here dest_pe = %d\n", dest_pe);
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

