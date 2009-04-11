TreePiece::TreePiece(nodeptr p, int which, bool isTopLevel_, int level_){
  isTopLevel = isTopLevel_;
  myLevel = level_;
  if(!isTopLevel_){
    // FIXME - myRoot should be set here.
    Parent(myRoot) = p;
    whichChildAmI = which;
  }
  numTotalMsgs = -1
  numRecvdMsgs = 0;
  haveChildren = false;
  myNumParticles = 0;

  haveCounts = false;

  for(int i = 0; i < NSUB; i++){
    sentTo[i] = 0;
  }
}

void TreePiece::recvTotalMsgCountsFromPieces(int totalNumFromParent){
  CkAssert(!isTopLevel);
  numTotalMsgs = totalNumFromParent;
  haveCounts = true;
  checkCompletion();
}


void TreePiece::recvTotalMsgCountsFromChunks(CkReductionMsg *msg){
  if(isTopLevel){
    int *data = (int *)msg->getData();
    numTotalMsgs = data[thisIndex]; /* between 0 and numTreePieces, 
                                       which is the number of 
                                       top-level treepieces */
    haveCounts = true;
    checkCompletion();
  }
  delete msg;
}

void TreePiece::checkCompletion(){
  if(numRecvdMsgs == numTotalMsgs){
    // the parent will not send any more messages
    if(haveChildren){
      // tell children that they will not receive any more messages 
      for(int i = 0; i < NSUB; i++){
        pieces[childrenTreePieces[i]].recvTotalMsgCountsFromPieces(sentTo[i]);
      }
    }
    else{
      // don't have children, build own tree
      buildTree();
    }
    // FIXME - contribute to reduction. need callback 
    // maincb to main
    contribute(0,0,CkReduction::concat,maincb);
  }
}

void TreePiece::recvParticles(ParticleMsg *msg){ 
  bodyptr particles = (bodyptr)msg->getData();
  numRecvdMsgs++;

  if(myNumParticles+msg->num > MAX_PARTS_PER_TP && !haveChildren){
    // insert children into pieces array
    for(int i = 0; i < NSUB; i++){
      int child = NSUB*thisIndex+numTreePieces+i;
      pieces[child].insert(myRoot, i, false, myLevel >> 1);
      childrenTreePieces[i] = child;
    }
    // FIXME - send own particles to children
    haveChildren = true;
    myNumParticles = 0;
  }

  if(haveChildren){
    CkVec<CkVec<bodyptr> > partsToChild;
    partsToChild.resize(NSUB);

    int num = msg->num;
    int xp[NDIM];
    bodyptr p;
    for(int i = 0; i < num; i++){
      int c; // part i goes to child c
      p = particles[i]; 
      c = NSUB*thisIndex+numTreePieces+subindex(xp,Pos(p));
      partsToChild[c].push_back(i);
    }

    // at this point, we have a list of particles 
    // destined for each child
    for(int c = 0; c < NSUB; c++){
      int len = partsToChild[c].length();
      if(len > 0){
        // create msg from partsToChild[c], send
        ParticleMsg *msg = new (len) ParticleMsg;
        msg->num = len;
        for(int i = 0; i < len; i++){
          bodyptr tmp = (partsToChild[c])[i];
          msg->particles[i] = *tmp;
        }
        sentTo[c]++;
        pieces[childrenTreePieces[c]].recvParticles(msg);
      }
    }
  }
  else{
    myNumParticles += msg->numParticles;
    // FIXME - add msg->particles to own particles
  }

  if(haveCounts){
    checkCompletion();
  }
  delete msg;
}

void TreePiece::buildTree(){
  // FIXME - build local tree here myParticles
  int l, xq[NDIM], xp[NDIM], xor[NDIM], subindex(), flag;
  int i, j, root_level;
  bool valid_root;
  int kidIndex;
  volatile nodeptr *volatile qptr, mynode;
  cellptr c;
  leafptr le;

  intcoord(xp, Pos(p));
  valid_root = TRUE;
  /*
  for (i = 0; i < NDIM; i++) {
    xor[i] = xp[i] ^ Local[ProcessId].Root_Coords[i];
  }
  for (i = IMAX >> 1; i > Level(root); i >>= 1) {
    for (j = 0; j < NDIM; j++) {
      if (xor[j] & i) {
        valid_root = FALSE;
        break;
      }
    }
    if (!valid_root) {
      break;
    }
  }
  if (!valid_root) {
    if (root != Global->G_root) {
      root_level = Level(root);
      for (j = i; j > root_level; j >>= 1) {
        root = (cellptr) Parent(root);
      }
      valid_root = TRUE;
      for (i = IMAX >> 1; i > Level(root); i >>= 1) {
        for (j = 0; j < NDIM; j++) {
          if (xor[j] & i) {
            valid_root = FALSE;
            break;
          }
        }
        if (!valid_root) {
          printf("P%d body %d\n", ProcessId, p - bodytab);
          root = Global->G_root;
        }
      }
    }
  }
  */
  root = Global->G_root;
  mynode = (nodeptr) root;
  kidIndex = subindex(xp, Level(mynode));
  qptr = &Subp(mynode)[kidIndex];

  l = Level(mynode) >> 1;

  flag = TRUE;
  while (flag) {                           /* loop descending tree     */
    if (l == 0) {
      error("not enough levels in tree\n");
    }
    if (*qptr == NULL) { 
      /* lock the parent cell */
      ALOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
      if (*qptr == NULL) {
        le = InitLeaf((cellptr) mynode, ProcessId);
        Parent(p) = (nodeptr) le;
        Level(p) = l;
        ChildNum(p) = le->num_bodies;
        ChildNum(le) = kidIndex;
        Bodyp(le)[le->num_bodies++] = p;
        *qptr = (nodeptr) le;
        flag = FALSE;
      }
      AULOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
      /* unlock the parent cell */
    }
    if (flag && *qptr && (Type(*qptr) == LEAF)) {
      /*   reached a "leaf"?      */
      ALOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
      /* lock the parent cell */
      if (Type(*qptr) == LEAF) {             /* still a "leaf"?      */
        le = (leafptr) *qptr;
        if (le->num_bodies == MAX_BODIES_PER_LEAF) {
          *qptr = (nodeptr) SubdivideLeaf(le, (cellptr) mynode, l,
              ProcessId);
        }
        else {
          Parent(p) = (nodeptr) le;
          Level(p) = l;
          ChildNum(p) = le->num_bodies;
          Bodyp(le)[le->num_bodies++] = p;
          flag = FALSE;
        }
      }
      AULOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
      /* unlock the node           */
    }
    if (flag) {
      mynode = *qptr;
      kidIndex = subindex(xp, l);
      qptr = &Subp(*qptr)[kidIndex];  /* move down one level  */
      l = l >> 1;                            /* and test next bit    */
    }
  }
  SETV(Local[ProcessId].Root_Coords, xp);
  return Parent((leafptr) *qptr);

}

