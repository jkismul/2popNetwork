COMMENT
From ModelDB entry 144401 by Nikita Vladimirov (https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=144401&file=/VladimirovTuTraub2012/gap.mod)
ENDCOMMENT


NEURON {
   POINT_PROCESS Gap
   POINTER vgap
   RANGE g, i   
   NONSPECIFIC_CURRENT i
}

PARAMETER { g = 1.0 (microsiemens) }

ASSIGNED {
   v (millivolt)
   vgap (millivolt)
   i (nanoamp)
}

BREAKPOINT { i = (v - vgap)*g }