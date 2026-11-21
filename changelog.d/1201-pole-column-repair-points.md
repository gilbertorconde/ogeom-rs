- **Fixed: slivers near a sphere's pole wound against the surface** (#1201).
  - A pole or an apex is one chart point at whatever `u` its ring gave it.
    The mesher's repair pass measured a triangle reaching it along a spiral
    from that `u`, and put the triangle's centre on another meridian. The
    triangles drawn round those centres folded over in space. The pass now
    measures and places along the column the triangle meets the pole by.
  - On `dome` and `dome_nurbs` at chord 0.01/3, 4 and 2 triangles faced
    inwards; now none do. Faces that reach a pole or an apex mesh in 1 to 9%
    fewer triangles, still within their chord. A face with no pole or apex
    meshes as it did.
  - Over the step-conformance parts at chords 0.01/3, 0.1 and 1.0, every
    triangle of a face that meets its chord is now wound with the face.
