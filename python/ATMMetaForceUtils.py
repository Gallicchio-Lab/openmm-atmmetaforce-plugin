%pythoncode %{

class ATMMetaForceUtils(object):

    def __init__(self, system, fix_zero_LJparams = True):
        self.system = system
        #place all existing forces in group 1, except for the non-bonded forces that go in group 2
        import re
        nbpattern = re.compile(".*Nonbonded.*")
        for force in system.getForces():
            if nbpattern.match(str(type(force))):
                force.setForceGroup(2)
                if fix_zero_LJparams:
                    self.fixZeroLJParams(force)
            else:
                force.setForceGroup(1)

    def _setRcptReferenceParticles(self, rcpt_ref_particles):
        if len(rcpt_ref_particles) != 3:
            raise ValueError("Invalid number of reference particles")
        self.rcpt_ref_particles = rcpt_ref_particles

    def _setLigReferenceParticles(self, lig_ref_particles):
        if len(lig_ref_particles) != 3:
            raise ValueError("Invalid number of reference particles")
        if not all(x in self.ligparticles for x in lig_ref_particles):
            raise ValueError("Invalid ligand atom indexes")
        self.lig_ref_particles = lig_ref_particles

    def addRestraintForce(self,
                          lig_cm_particles = None, rcpt_cm_particles = None,
                          kfcm = 0.0 * kilocalorie_per_mole/angstrom**2,
                          tolcm = 0.0 * angstrom,
                          lig_ref_particles = None, rcpt_ref_particles = None,
                          angle_center = 0.*degrees,
                          kfangle = 0.0 * kilocalorie_per_mole/degrees**2,
                          angletol = 10.0*degrees,
                          dihedral1center = 0.*degrees,
                          kfdihedral1 = 0.0 * kilocalorie_per_mole/degrees**2,
                          dihedral1tol = 10.0*degrees,
                          dihedral2center = 0.*degrees,
                          kfdihedral2 = 0.0 * kilocalorie_per_mole/degrees**2,
                          dihedral2tol = 10.0*degrees,
                          offset = [0., 0., 0.] * angstrom):

        if not (lig_cm_particles and rcpt_cm_particles):
            return None

        if not (len(lig_cm_particles) > 0 and len(rcpt_cm_particles) > 0):
            return None

        do_angles = False
        if lig_ref_particles and rcpt_ref_particles:
            if not (len(lig_ref_particles) == 3 and len(rcpt_ref_particles) == 3):
                raise ValueError("Invalid lists of reference atoms")
            do_angles = True

        expr = ""

        expr += "1.0 *( (kfcm/2)*step(d12-tolcm)*(d12-tolcm)^2 )                      "

        if do_angles:
            expr += " + 1.0 *( (kfcd0/2)*(step(dm0)*max(0,db0)^2+step(-dm0)*max(0,-da0)^2) )  "
            expr += " + 1.0 *( (kfcd1/2)*(step(dm1)*max(0,db1)^2+step(-dm1)*max(0,-da1)^2) )  "
            expr += " + 1.0 *( (kfcd2/2)*(step(dm2)*max(0,db2)^2+step(-dm2)*max(0,-da2)^2) )  "

        expr += " ; d12 = sqrt((x1 - offx - x2)^2 + (y1 - offy - y2)^2 + (z1 - offz - z2)^2 ) ; "

        if do_angles:
            expr += "db0 = xb0 - pi*floor(xb0/pi + 0.5)  ; xb0 = theta - b0 ; "
            expr += "da0 = xa0 - pi*floor(xa0/pi + 0.5)  ; xa0 = theta - a0 ; "
            expr += "dm0 = xm0 - pi*floor(xm0/pi + 0.5)  ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; theta = angle(g3,g6,g7) ; "

            expr += "db1 = xb1 - twopi*floor(xb1/twopi + 0.5)  ; xb1 = phi1 - b1 ; "
            expr += "da1 = xa1 - twopi*floor(xa1/twopi + 0.5)  ; xa1 = phi1 - a1 ; "
            expr += "dm1 = xm1 - twopi*floor(xm1/twopi + 0.5)  ; xm1 = phi1 - mid1 ; mid1 = (a1 + b1)/2 ; phi1 = dihedral(g4,g3,g6,g7) ; "

            expr += "db2 = xb2 - twopi*floor(xb2/twopi + 0.5)  ; xb2 = phi2 - b2 ; "
            expr += "da2 = xa2 - twopi*floor(xa2/twopi + 0.5)  ; xa2 = phi2 - a2 ; "
            expr += "dm2 = xm2 - twopi*floor(xm2/twopi + 0.5)  ; xm2 = phi2 - mid2 ; mid2 = (a2 + b2)/2 ; phi2 = dihedral(g3,g6,g7,g8) ; "

        expr += "twopi = 2*pi ;"
        expr += "pi = %f" % math.pi

        if do_angles:
            force =  mm.CustomCentroidBondForce(8,expr)
        else:
            force =  mm.CustomCentroidBondForce(2,expr)

        self.system.addForce(force)

        force.setForceGroup(1) #the restraint force will be evaluated separately

        force.addPerBondParameter("kfcm")
        force.addPerBondParameter("tolcm")
        force.addPerBondParameter("offx")
        force.addPerBondParameter("offy")
        force.addPerBondParameter("offz")

        force.addGroup(lig_cm_particles) #g1 CM of lig
        force.addGroup(rcpt_cm_particles) #g2 CM of rcpt

        if do_angles:

            force.addPerBondParameter("kfcd0")
            force.addPerBondParameter("a0")
            force.addPerBondParameter("b0")

            force.addPerBondParameter("kfcd1")
            force.addPerBondParameter("a1")
            force.addPerBondParameter("b1")

            force.addPerBondParameter("kfcd2")
            force.addPerBondParameter("a2")
            force.addPerBondParameter("b2")

            force.addGroup([rcpt_ref_particles[0]]) #g3 rcpt ref 0
            force.addGroup([rcpt_ref_particles[1]]) #g4 rcpt ref 1
            force.addGroup([rcpt_ref_particles[2]]) #g5 rcpt ref 2

            force.addGroup([lig_ref_particles[0]]) #g6 lig ref 0
            force.addGroup([lig_ref_particles[1]]) #g7 lig ref 1
            force.addGroup([lig_ref_particles[2]]) #g8 lig ref 2

        kfc = kfcm / (kilojoule_per_mole/radians**2)
        tolc = tolcm / nanometer
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]

        if do_angles:

            kfcd0 = kfangle / (kilojoule_per_mole/radians**2)
            a0 = (angle_center - angletol)/radians
            b0 = (angle_center + angletol)/radians

            kfcd1 = kfdihedral1 / (kilojoule_per_mole/radians**2)
            a1 = (dihedral1center - dihedral1tol)/radians
            b1 = (dihedral1center + dihedral1tol)/radians

            kfcd2 = kfdihedral2 / (kilojoule_per_mole/radians**2)
            a2 = (dihedral2center - dihedral2tol)/radians
            b2 = (dihedral2center + dihedral2tol)/radians

            groups = range(8)
            params = [kfc, tolc, offx, offy, offz, kfcd0, a0, b0, kfcd1, a1, b1, kfcd2, a2, b2]
        else:
            groups = [0,1]
            params = [kfc, tolc, offx, offy, offz]

        force.addBond(groups, params)
        return force

    # a force to keep two ligands aligned
    # keeps ref atoms 1 near each other when offset is added
    #    lig2 is assumed to be translated by offset relative to lig1
    #     (lig2pos = lig1pos + offset)
    # keeps theta angle near zero
    # keeps psi angle near zero (assuming theta is near zero)
    def addAlignmentForce(self,
                          liga_ref_particles = None, ligb_ref_particles = None,
                          kfdispl = 0.0 * kilocalorie_per_mole/angstrom**2,
                          ktheta = 0.0 * kilocalorie_per_mole,
                          kpsi = 0.0 * kilocalorie_per_mole,
                          offset = [0., 0., 0.] * angstrom):

        if liga_ref_particles and ligb_ref_particles:
            if not (len(liga_ref_particles) == 3 and len(ligb_ref_particles) == 3):
                raise ValueError("Invalid lists of reference atoms")

        expr  = " (kfdispl/2)*distsq  ; " #displacement potential
        expr += " distsq = (x1 - offx - x2)^2 + "
        expr += "          (y1 - offy - y2)^2 + "
        expr += "          (z1 - offz - z2)^2   " #square distance between b1 and a1 after displacing back b

        displforce =  mm.CustomCompoundBondForce(2, expr);
        self.system.addForce(displforce)
        displforce.addPerBondParameter("kfdispl")
        displforce.addPerBondParameter("offx")
        displforce.addPerBondParameter("offy")
        displforce.addPerBondParameter("offz")
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]
        displforce.addBond([ligb_ref_particles[0],liga_ref_particles[0]] ,
                           [kfdispl/((kilojoule_per_mole/nanometer**2)),
                            offx, offy, offz])
        displforce.setForceGroup(1);


        expr = "(ktheta/2) * (1 - cost) ;" #theta restraint

        expr += "cost = xdn1*xdn2 + ydn1*ydn2 + zdn1*zdn2 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 ;"
        expr += "xdn2 = xd2/dn2 ; ydn2 = yd2/dn2 ; zdn2 = zd2/dn2 ;"
        expr += "dn2 = sqrt(xd2^2 + yd2^2 + zd2^2 ) ;"
        expr += "xd2 = x4 - x3 ; "
        expr += "yd2 = y4 - y3 ; "
        expr += "zd2 = z4 - z3   "

        thetaforce = mm.CustomCompoundBondForce(4, expr)
        self.system.addForce(thetaforce)
        thetaforce.addPerBondParameter("ktheta")
        thetaforce.addBond([ligb_ref_particles[0],ligb_ref_particles[1],liga_ref_particles[0],liga_ref_particles[1] ] ,
                           [ktheta/kilojoule_per_mole])
        thetaforce.setForceGroup(1);

        expr = "(kpsi/2) * (1 - cosp) ;" #psi restraint

        expr += "cosp = xvn*xwn + yvn*ywn + zvn*zwn ; "
        expr += "xvn = xv/v ; yvn = yv/v; zvn = zv/v ;"
        expr += "v = sqrt(xv^2 + yv^2 + zv^2 ) ;"
        expr += "xv = xd0 - dot01*xdn1 ;"
        expr += "yv = yd0 - dot01*ydn1 ;"
        expr += "zv = zd0 - dot01*zdn1 ;"
        expr += "dot01 =  xd0*xdn1 +  yd0*ydn1 +  zd0*zdn1 ;"
        expr += "xd0 = x3 - x1 ;"
        expr += "yd0 = y3 - y1 ;"
        expr += "zd0 = z3 - z1 ;"
        expr += "xwn = xw/w ; ywn = yw/w; zwn = zw/w ;"
        expr += "w = sqrt(xw^2 + yw^2 + zw^2 ) ;"
        expr += "xw = xd3 - dot31*xdn1 ;"
        expr += "yw = yd3 - dot31*ydn1 ;"
        expr += "zw = zd3 - dot31*zdn1 ;"
        expr += "dot31 =  xd3*xdn1 +  yd3*ydn1 +  zd3*zdn1 ;"
        expr += "xd3 = x5 - x4 ;"
        expr += "yd3 = y5 - y4 ;"
        expr += "zd3 = z5 - z4 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 "

        psiforce = mm.CustomCompoundBondForce(5, expr)
        self.system.addForce(psiforce)
        psiforce.addPerBondParameter("kpsi")
        psiforce.addBond([ligb_ref_particles[0],ligb_ref_particles[1],ligb_ref_particles[2],
                          liga_ref_particles[0],liga_ref_particles[2] ] ,
                           [kpsi/kilojoule_per_mole])
        psiforce.setForceGroup(1)
        return (displforce, thetaforce, psiforce )

    # applies flat-bottom restraints to a set of atoms based on a set of reference positions
    def addPosRestraints(self, particles, refpos, fc = 25.0 * kilocalorie_per_mole/angstrom**2, tol = 0.5 * angstrom, periodic = True):
        if not particles or len(particles) == 0:
            return None
        if periodic:
            posrestforce = mm.CustomExternalForce("0.5*fc*select(step(dist-tol), (dist-tol)^2, 0); dist = periodicdistance(x,y,z,x0,y0,z0)")
        else:
            posrestforce = mm.CustomExternalForce("0.5*fc*select(step(dist-tol), (dist-tol)^2, 0); dist = sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)")

        posrestforce.addPerParticleParameter("x0")
        posrestforce.addPerParticleParameter("y0")
        posrestforce.addPerParticleParameter("z0")
        posrestforce.addPerParticleParameter("fc")
        posrestforce.addPerParticleParameter("tol")

        posrestforce.setForceGroup(1)

        self.system.addForce(posrestforce)

        for p in particles:
            x0 = refpos[p][0]
            y0 = refpos[p][1]
            z0 = refpos[p][2]
            posrestforce.addParticle(p, [x0, y0, z0, fc, tol])
        return posrestforce

    #fix zero LJs
    def fixZeroLJParams(self, force):
        small = 1.0e-6
        if isinstance(force, mm.NonbondedForce):
            for i in range(force.getNumParticles()):
                nbparam = force.getParticleParameters(i)
                charge = nbparam[0]
                sigmaLJ = nbparam[1]
                epsilonLJ = nbparam[2]
                if sigmaLJ._value < small and epsilonLJ._value < small:
                    sigmaLJ = 0.1 * angstrom
                    epsilonLJ = 1.e-4 * kilocalories_per_mole
                    force.setParticleParameters(i, charge, sigmaLJ, epsilonLJ)

%}
