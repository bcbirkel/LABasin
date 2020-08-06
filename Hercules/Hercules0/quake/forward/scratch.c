// To debug 
  for ( iPlane = 0; iPlane < theNumberOfPlanes; iPlane++ ){
    fprintf(stdout,"\n MyID= %d SAS=%d SDD=%d TI=%d",myID,thePlanes[iPlane].numberofstepsalongstrike
 	                                                 ,thePlanes[iPlane].numberofstepsdowndip
                                                         ,thePlanes[iPlane].typeplaneinput );
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

 MPI_Abort(MPI_COMM_WORLD, ERROR);



 // to debug output the broadcast info to check that is the same read per pe
  if(myID==3){
    FILE *fpT;
    fpT=fopen("testInput_pe3","w");

    iPlane=0;
    nI=thePlanes[iPlane].numberofstepsalongstrike;
    nJ=thePlanes[iPlane].numberofstepsdowndip;
    for (iC=0; iC < (nI*nJ); iC++)   
      fprintf(fpT," \n %lf %lf %lf ",thePlanes[iPlane].latIn[iC],thePlanes[iPlane].lonIn[iC],
	      thePlanes[iPlane].depthIn[iC]);
    
    fclose(fpT);
  }









	//}
	//    }
	//}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	MPI_Abort(MPI_COMM_WORLD, ERROR);
	//	for ( iStrike = 0; iStrike < thePlanes[iPlane].numberofstepsalongstrike; iStrike++ ){
	//xLocal = iStrike*thePlanes[iPlane].stepalongstrike;
	//	for ( iDownDip = 0; iDownDip < thePlanes[iPlane].numberofstepsdowndip; iDownDip++ ) {








