- **Fixed: a round whose ball is tangent to a far face along its whole
  run no longer hangs** (#1171). The far-edge clearance proof halved each
  interval to depth 24 on both sides when the face stood exactly a radius
  away, about 16 million distance measures per station. It now stops after
  2^17 measures per interval, about three times the most any passing proof
  uses, and treats the face as not shown clear. A slot exactly twice the
  radius wide used to hang. Now its two floor edges give the half-circle
  trough, and one edge alone gives its band.
