__kernel void HybridForce(int numParticles,
			  __global real4* restrict force,
			  __global real4* restrict force_state1,
			  __global real4* restrict force_state2,
			  float sp
			  ){
  uint i = get_global_id(0);
  real lmb  = sp;
  real lmb1 = 1.0f - sp;
  while (i < numParticles) {
    force[i] += lmb*force_state2[i]+ lmb1*force_state1[i];
    i += get_global_size(0);
  }
}

__kernel void CopyState(int numParticles,
			__global real4* restrict posq,
			__global real4* restrict posq1,
			__global real4* restrict posq2,
			__global float4* restrict displ
#ifdef USE_MIXED_PRECISION
			,
			__global real4* restrict posqCorrection,
			__global real4* restrict posq1Correction,
			__global real4* restrict posq2Correction
#endif
			){

  //set the coordinates of the context for state 1
  int i = get_global_id(0);
  while (i < numParticles) {
    posq1[i] = posq[i];
#ifdef USE_MIXED_PRECISION
    posq1Correction[i] = posqCorrection[i];
#endif
    i += get_global_size(0);
  }
  
  //set the coordinates of the context for state 2
  i = get_global_id(0);
  while (i < numParticles) {
    real4 d = (real4)(displ[i].x, displ[i].y, displ[i].z, 0);
    posq2[i] = posq[i] + d;
#ifdef USE_MIXED_PRECISION
    posq2Correction[i] = posqCorrection[i];
#endif
    i += get_global_size(0);
  }
}


