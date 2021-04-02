__kernel void HybridForce(int numParticles,
			  __global real4* restrict force,
			  __global real4* restrict force_state1,
			  __global real4* restrict force_state2,
			  float sp,
			  __global int* restrict atomOrder,
			  __global int* restrict invAtomOrder1,
			  __global int* restrict invAtomOrder2
			  ){
  uint i = get_global_id(0);
  real lmb  = sp;
  real lmb1 = 1.0f - sp;
  while (i < numParticles) {//here i is the index in the outer context
    int i1 = invAtomOrder1[atomOrder[i]];
    int i2 = invAtomOrder2[atomOrder[i]];
    force[i] += lmb*force_state2[i2]+ lmb1*force_state1[i1];
    i += get_global_size(0);
  }
}

__kernel void CopyState(int numParticles,
			__global real4* restrict posq,
			__global real4* restrict posq1,
			__global real4* restrict posq2,
			__global float4* restrict displ,
			__global int* restrict atomOrder,
			__global int* restrict invAtomOrder1,
			__global int* restrict invAtomOrder2
			){

  //set the coordinates of the context for state 1
  int i = get_global_id(0);
  while (i < numParticles) { //here i is the index in the outer context
    int i1 = invAtomOrder1[atomOrder[i]];
    posq1[i1] = posq[i];
    i += get_global_size(0);
  }
  
  //set the coordinates of the context for state 2
  i = get_global_id(0);
  while (i < numParticles) { //outer context's index
    int index = atomOrder[i]; //system's index
    int i2 = invAtomOrder2[index]; //inner context 2 index
    real4 d = (real4)(displ[index].x, displ[index].y, displ[index].z, 0);
    posq2[i2] = posq[i] + d;
    i += get_global_size(0);
  }
}


