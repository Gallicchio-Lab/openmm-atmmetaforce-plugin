__kernel void HybridForce(int numParticles,
			  __global real4* restrict force,
			  __global real4* restrict force_state1,
			  __global real4* restrict force_state2,
			 float sp){
  uint index = get_global_id(0);
  real lmb  = sp;
  real lmb1 = 1.0f - sp;
  while (index < numParticles) {
    force[index] += lmb*force_state2[index]+ lmb1*force_state1[index];
    index += get_global_size(0);
  }
}

__kernel void CopyState(int numParticles,
			  __global real4* restrict posq,
			  __global real4* restrict posq1,
			  __global real4* restrict posq2,
			  __global float4* restrict displ){
  uint index = get_global_id(0);
  while (index < numParticles) {
    posq1[index] = posq[index];
    index += get_global_size(0);
  }
  index = get_global_id(0);
  while (index < numParticles) {
    posq2[index] = posq[index];
    index += get_global_size(0);
  }
  index = get_global_id(0);
  while (index < numParticles) {
    posq2[index] += (real4)(displ[index].xyz,0);
    index += get_global_size(0);
  }
}


