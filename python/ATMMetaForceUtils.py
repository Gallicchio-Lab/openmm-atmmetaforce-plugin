%pythoncode %{
import string
import random
import numpy as np
from simtk.openmm import version as ommversion

class ATMMetaForceUtils(object):

    def __init__(self, system, fix_zero_LJparams = True):
        self.system = system

        self.CMCMDistForce = None
        self.CMAngleThetaForce = None
        self.CMAnglePhiForce = None
        self.CMAnglePsiForce = None

        self.major_ommversion = int(ommversion.version.split(".")[0])
        self.minor_ommversion = int(ommversion.version.split(".")[1])

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
                          offset = [0., 0., 0.] * angstrom):
        print("warning: AddRestraintForce() is deprecated. Use addVsiteRestraintForceCMCM()")
        res = self.addVsiteRestraintForceCMCM(lig_cm_particles, rcpt_cm_particles,
                                              kfcm, tolcm, offset)
        return res


    #
    # Adds restrains to keep the the CM of the ligand atoms around the CM of the receptor atoms
    # plus an offset, within a tolerance, with a flat-bottom harmonic potential
    def addVsiteRestraintForceCMCM(self,
                          lig_cm_particles = None, rcpt_cm_particles = None,
                          kfcm = 0.0 * kilocalorie_per_mole/angstrom**2,
                          tolcm = 0.0 * angstrom,
                          offset = [0., 0., 0.] * angstrom):
        assert len(lig_cm_particles) > 0
        assert len(rcpt_cm_particles) > 0

        expr = "(kfcm/2)*step(d12-tolcm)*(d12-tolcm)^2 "
        expr += " ; d12 = sqrt((x1 - offx - x2)^2 + (y1 - offy - y2)^2 + (z1 - offz - z2)^2 ) ; "
        force = None
        if self.CMCMDistForce is None:
           self.CMCMDistForce  =  mm.CustomCentroidBondForce(2,expr)
           force = self.CMCMDistForce
           force.addPerBondParameter("kfcm")
           force.addPerBondParameter("tolcm")
           force.addPerBondParameter("offx")
           force.addPerBondParameter("offy")
           force.addPerBondParameter("offz")
           numgroups = 0
           force.setForceGroup(1)
           self.system.addForce(force)
        else:
            force = self.CMCMDistForce
            numgroups = force.getNumGroups()
        force.addGroup(lig_cm_particles) #g1 CM of lig
        force.addGroup(rcpt_cm_particles) #g2 CM of rcpt
        kfc = kfcm / (kilojoule_per_mole/nanometer**2)
        tolc = tolcm / nanometer
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]
        groups = [numgroups + 0, numgroups + 1]
        params = [kfc, tolc, offx, offy, offz]
        force.addBond(groups, np.array(params, dtype=np.double))
        return force

    def _roundExpression(self, n, x):
        return "{n} = select(step({x} - floor({x})  - 0.5), ceil({x}), floor({x}))".format(n = n, x = x)

    def _wrapExpression(self, w, x, r):
        letters = string.ascii_lowercase
        n = "n" + ''.join(random.choice(letters) for i in range(4))
        expr = "{w} = {x} - {r}*{n} ; ".format(w = w, x = x, r = r, n = n) + " ; "
        expr += self._roundExpression("{n}".format(n = n), "({x}/{r})".format(x = x, r = r))
        return expr

    def _dotExpression(self, dot, x1, y1, z1, x2, y2, z2):
        return "{dot} = {x1}*{x2} + {y1}*{y2} + {z1}*{z2}".format(dot = dot, x1 = x1, y1 = y1, z1 = z1,  x2 = x2, y2 = y2, z2 = z2) 

    def _unitvExpression(self, xu, yu, zu, x, y, z):
        letters = string.ascii_lowercase
        d = "d" + ''.join(random.choice(letters) for i in range(4))
        return "{xu} = {x}/{d} ; {yu} = {y}/{d} ; {zu} = {z}/{d} ; {d} = sqrt({x}*{x}+{y}*{y}+{z}*{z})".format(
            d = d, xu = xu, yu = yu, zu = zu, x = x, y = y, z = z)

    def _diffvExpression(self, dx, dy, dz, x2, y2, z2, x1, y1, z1):
        return "{dx} = {x2} - {x1} ; {dy} = {y2} - {y1} ; {dz} = {z2} - {z1} ".format(
            dx = dx, dy = dy, dz = dz,
            x2 = x2, y2 = y2, z2 = z2,
            x1 = x1, y1 = y1, z1 = z1)

    def _crossExpression(self, xc, yc, zc, x1, y1, z1, x2, y2, z2):
        return ( "{xc} =  {y1}*{z2} - {y2}*{z1} ; "
                 "{yc} = -{x1}*{z2} + {x2}*{z1} ; "
                 "{zc} =  {x1}*{y2} - {x2}*{y1}   " ).format(xc = xc, yc = yc, zc = zc, x1 = x1, y1 = y1, z1 = z1, x2 = x2, y2 = y2, z2 = z2)

    def _cosangleExpression(self, costh, x1, y1, z1, x2, y2, z2, x3, y3, z3):
        if self.major_ommversion >= 7 and self.minor_ommversion >= 6:
            return costh + "= cos(pointangle(" + x1 + "," + y1 + "," + z1 + "," + \
                                                 x2 + "," + y2 + "," + z2 + "," + \
                                                 x3 + "," + y3 + "," + z3 + "))"
        else:
            letters = string.ascii_lowercase
            v1  = "v1"  + ''.join(random.choice(letters) for i in range(4))
            v1u = "v1u" + ''.join(random.choice(letters) for i in range(4))
            v2  = "v2"  + ''.join(random.choice(letters) for i in range(4))
            v2u = "v2u" + ''.join(random.choice(letters) for i in range(4))
            expr  = self._dotExpression(costh, "x"+v1u,"y"+v1u,"z"+v1u,"x"+v2u,"y"+v2u,"z"+v2u) + " ; "
            expr += self._unitvExpression("x"+v1u,"y"+v1u,"z"+v1u, "x"+v1,"y"+v1,"z"+v1) + " ; "
            expr += self._unitvExpression("x"+v2u,"y"+v2u,"z"+v2u, "x"+v2,"y"+v2,"z"+v2) + " ; "
            expr += self._diffvExpression("x"+v1, "y"+v1, "z"+v1, x1, y1, z1, x2, y2, z2) + " ; "
            expr += self._diffvExpression("x"+v2, "y"+v2, "z"+v2, x3, y3, z3, x2, y2, z2)
            return expr

    def _dihedralExpression(self, phi, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
        if self.major_ommversion >= 7 and self.minor_ommversion >= 6:
            return phi + "= pointdihedral(" + x1 + "," + y1 + "," + z1 + "," + \
                                              x2 + "," + y2 + "," + z2 + "," + \
                                              x3 + "," + y3 + "," + z3 + "," + \
                                              x4 + "," + y4 + "," + z4 + ")"
        else:
            letters = string.ascii_lowercase
            v1   = "v1"   + ''.join(random.choice(letters) for i in range(4))
            v2   = "v2"   + ''.join(random.choice(letters) for i in range(4))
            v2u  = "v2u"  + ''.join(random.choice(letters) for i in range(4))
            v3   = "v3"   + ''.join(random.choice(letters) for i in range(4))
            n1   = "n1"   + ''.join(random.choice(letters) for i in range(4))
            n1u  = "n1u"  + ''.join(random.choice(letters) for i in range(4))
            n2   = "n2"   + ''.join(random.choice(letters) for i in range(4))
            n2u  = "n2u"  + ''.join(random.choice(letters) for i in range(4))
            m1   = "m1"   + ''.join(random.choice(letters) for i in range(4))
            m1u  = "m1u"  + ''.join(random.choice(letters) for i in range(4))
            cx   = "cx"   + ''.join(random.choice(letters) for i in range(4))
            cy   = "cy"   + ''.join(random.choice(letters) for i in range(4))
            expr  = "{phi} = atan2({cy} , {cx}) ; ".format(phi = phi, cx = cx, cy = cy)
            expr += self._dotExpression(cx, "x"+n1u,"y"+n1u,"z"+n1u,"x"+n2u,"y"+n2u,"z"+n2u) + " ; "
            expr += self._dotExpression(cy, "x"+m1u,"y"+m1u,"z"+m1u,"x"+n2u,"y"+n2u,"z"+n2u) + " ; "
            expr += self._unitvExpression("x"+m1u,"y"+m1u,"z"+m1u, "x"+m1,  "y"+m1,  "z"+m1) + " ; "
            expr += self._crossExpression("x"+m1, "y"+m1, "z"+m1,  "x"+v2u, "y"+v2u, "z"+v2u, "x"+n1u, "y"+n1u, "z"+n1u) + " ; "
            expr += self._unitvExpression("x"+n1u,"y"+n1u,"z"+n1u, "x"+n1,  "y"+n1,  "z"+n1) + " ; "
            expr += self._unitvExpression("x"+n2u,"y"+n2u,"z"+n2u, "x"+n2,  "y"+n2,  "z"+n2) + " ; "
            expr += self._unitvExpression("x"+v2u,"y"+v2u,"z"+v2u, "x"+v2,  "y"+v2,  "z"+v2) + " ; "
            expr += self._crossExpression("x"+n2, "y"+n2, "z"+n2,  "x"+v2,  "y"+v2,  "z"+v2, "x"+v3, "y"+v3, "z"+v3) + " ; "
            expr += self._crossExpression("x"+n1, "y"+n1, "z"+n1,  "x"+v1,  "y"+v1,  "z"+v1, "x"+v2, "y"+v2, "z"+v2) + " ; "
            expr += self._diffvExpression("x"+v1, "y"+v1, "z"+v1, x2, y2, z2, x1, y1, z1) + " ; "
            expr += self._diffvExpression("x"+v2, "y"+v2, "z"+v2, x3, y3, z3, x2, y2, z2) + " ; "
            expr += self._diffvExpression("x"+v3, "y"+v3, "z"+v3, x4, y4, z4, x3, y3, z3)
            return expr

    # ** UNTESTED **
    # Adds Boresch-style Vsite restraints [J. Phys. Chem. B, Vol. 107, No. 35, 2003]
    # between 3 reference atoms of the ligand and
    # 3 reference atoms of the receptor with flat-bottom potentials
    # Reference particles of the receptor: a, b, c
    # Reference particles of the ligand:   A, B, C
    # rA = distance(A, a)
    # thetaA = angle(b, a, A)
    # phiA = dihedral(c, b, a, A)
    # (rA, thetaA, phiA are the spherical polar coordinates of A wrt the reference frame of the recepotor)
    # thetaB = angle(a, A, B)
    # phiB = dihedral(b, a, A, B)
    # phiC = dihedral(a, A, B, C)
    def _addVsiteRestraintForceBoresch(self,
                                      lig_ref_particles, rcpt_ref_particles,
                                      kfrA,  rA0,  rAtol,
                                      kfthA, thA0, thAtol,
                                      kfphA, phA0, phAtol,
                                      kfthB, thB0, thBtol,
                                      kfphB, phB0, phBtol,
                                      kfphC, phC0, phCtol):
        assert len( lig_ref_particles) == 3
        assert len(rcpt_ref_particles) == 3

        bondforce = None
        if kfrA is not None:
            bondforce = mm.CustomBondForce("0.5*kf*( step(d)*max(0,d-tol)^2 + step(-d)*max(0,-d-tol)^2 ) ; d = r - r0")
            bondforce.addPerBondParameter("kf")
            bondforce.addPerBondParameter("r0")
            bondforce.addPerBondParameter("tol")
            bondforce.setForceGroup(1)
            p0 = rcpt_ref_particles[0]
            p1 =  lig_ref_particles[0]
            kf = kfrA/(kilojoule_per_mole/nanometer**2)
            r0 = rA0/nanometer
            tol = rAtol/nanometer
            bondforce.addBond([p0,p1], np.array([kf, r0, tol ], dtype=np.double))
            self.system.addForce(bondforce)

        angleforce = None
        if kfthA is not None or kfthB is not None:
            expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
            expr += "db0 = xb0 - pi*floor(xb0/pi + 0.5)  ; xb0 = theta - b0 ; "
            expr += "da0 = xa0 - pi*floor(xa0/pi + 0.5)  ; xa0 = theta - a0 ; "
            expr += "dm0 = xm0 - pi*floor(xm0/pi + 0.5)  ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; "
            expr += "twopi = 2*pi ;"
            expr += "pi = %f" % math.pi
            angleforce = mm.CustomAngleForce(expr)
            angleforce.addPerAngleParameter("kf")
            angleforce.addPerAngleParameter("a0")
            angleforce.addPerAngleParameter("b0")
            angleforce.setForceGroup(1)
            if kfthA is not None:
                kf = kfthA/(kilojoule_per_mole/radians**2)
                a0 = (thA0 - thAtol)/radians
                b0 = (thA0 + thAtol)/radians
                p0 = rcpt_ref_particles[1]
                p1 = rcpt_ref_particles[0]
                p2 =  lig_ref_particles[0]
                angleforce.addAngle([p0,p1,p2], np.array([kf,a0,b0], dtype=np.double) )
            if kfthB is not None:
                kf = kfthB/(kilojoule_per_mole/radians**2)
                a0 = (thB0 - thBtol)/radians
                b0 = (thB0 + thBtol)/radians
                p0 = rcpt_ref_particles[0]
                p1 =  lig_ref_particles[0]
                p2 =  lig_ref_particles[1]
                angleforce.addAngle([p0,p1,p2], np.array([kf,a0,b0], dtype=np.double) )
            self.system.addForce(angleforce)

        torsforce = None
        if kfphA is not None or kfphB is not None or kfphC is not None:
            expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
            expr += "db0 = xb0 - twopi*floor(xb0/twopi + 0.5)  ; xb0 = theta - b0 ; "
            expr += "da0 = xa0 - twopi*floor(xa0/twopi + 0.5)  ; xa0 = theta - a0 ; "
            expr += "dm0 = xm0 - twopi*floor(xm0/twopi + 0.5)  ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; "
            expr += "twopi = 2*pi ;"
            expr += "pi = %f" % math.pi
            torsforce = mm.CustomTorsionForce(expr)
            torsforce.addPerTorsionParameter("kf")
            torsforce.addPerTorsionParameter("a0")
            torsforce.addPerTorsionParameter("b0")
            torsforce.setForceGroup(1)
            if kfphA is not None:
                kf = kfphA/(kilojoule_per_mole/radians**2)
                a0 = (phA0 - phAtol)/radians
                b0 = (phA0 + phAtol)/radians
                p0 = rcpt_ref_particles[2]
                p1 = rcpt_ref_particles[1]
                p2 = rcpt_ref_particles[0]
                p3 =  lig_ref_particles[0]
                torsforce.addTorsion([p0,p1,p2,p3], np.array([kf,a0,b0], dtype=np.double) )
            if kfphB is not None:
                kf = kfphB/(kilojoule_per_mole/radians**2)
                a0 = (phB0 - phBtol)/radians
                b0 = (phB0 + phBtol)/radians
                p0 = rcpt_ref_particles[1]
                p1 = rcpt_ref_particles[0]
                p2 =  lig_ref_particles[0]
                p3 =  lig_ref_particles[1]
                torsforce.addTorsion([p0,p1,p2,p3], np.array([kf,a0,b0], dtype=np.double) )
            if kfphC is not None:
                kf = kfphC/(kilojoule_per_mole/radians**2)
                a0 = (phC0 - phCtol)/radians
                b0 = (phC0 + phCtol)/radians
                p0 = rcpt_ref_particles[0]
                p1 =  lig_ref_particles[0]
                p2 =  lig_ref_particles[1]
                p3 =  lig_ref_particles[2]
                torsforce.addTorsion([p0,p1,p2,p3], np.array([kf,a0,b0], dtype=np.double) )
            self.system.addForce(torsforce)

        return (bondforce, angleforce, torsforce)

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
                           np.array([kfdispl/((kilojoule_per_mole/nanometer**2)),
                                     offx, offy, offz], dtype=np.double) )
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
                           np.array([ktheta/kilojoule_per_mole], dtype=np.double) )
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
                         np.array([0.5*kpsi/kilojoule_per_mole], dtype=np.double) )
        psiforce.addBond([liga_ref_particles[0],liga_ref_particles[1],liga_ref_particles[2],
                          ligb_ref_particles[0],ligb_ref_particles[2] ] ,
                         np.array([0.5*kpsi/kilojoule_per_mole], dtype=np.double) )
        psiforce.setForceGroup(1)
        return (displforce, thetaforce, psiforce )

    # restrain angle between alignment z-axes defined by centroids
    # similar to AlignmentForce but for centroids
    def addVsiteRestraintForceCMAngles(self,
                    lig_cm_groups = None, rcpt_cm_groups = None,
                    ktheta = None, theta0 = None, thetatol = None,
                    kphi = None, phi0 = None, phitol = None, 
                    kpsi = None, psi0 = None, psitol = None):
        assert len(lig_cm_groups) == 3
        assert len(rcpt_cm_groups) == 3

        #theta restraint
        # (kf/2) (cos(theta) - cos(theta0))^2 with a tolerance
        thetaforce = None
        expr = "0.5*kf*( step(dcost)*max(0,dcost-ctol)^2 + step(-dcost)*max(0,-dcost-ctol)^2 ) ; dcost = cost - cos0 ; "        
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
        if ktheta is not None:
            if self.CMAngleThetaForce is None:
               self.CMAngleThetaForce = mm.CustomCentroidBondForce(4, expr)
               thetaforce = self.CMAngleThetaForce
               numthetagroups = 0
               thetaforce.addPerBondParameter("kf")
               thetaforce.addPerBondParameter("cos0")
               thetaforce.addPerBondParameter("ctol")
               thetaforce.setForceGroup(1)
               self.system.addForce(thetaforce)
            else:
                thetaforce = self.CMAngleThetaForce
                numthetagroups = thetaforce.getNumGroups()
            thetaforce.addGroup(rcpt_cm_groups[0])
            thetaforce.addGroup(rcpt_cm_groups[1])
            thetaforce.addGroup( lig_cm_groups[0])
            thetaforce.addGroup( lig_cm_groups[1])
            kf = ktheta/(kilojoule_per_mole)
            cos0 = math.cos(theta0/radians)
            a0 = max(0, (theta0-thetatol)/radians)
            b0 = min(math.pi, (theta0+thetatol)/radians)
            ctol = abs(math.cos(a0) - math.cos(b0))
            groupids = [ numthetagroups + i for i in range(4) ]
            thetaforce.addBond(groupids, np.array([kf, cos0, ctol], dtype=np.double )  )

        #phi restraint
        # (kf/2) (phi - phi0)^2 with a tolerance and periodicity
        phiforce = None
        expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
        expr += self._wrapExpression("db0", "xb0", "twopi") + " ; xb0 = psi - b0 ; "
        expr += self._wrapExpression("da0", "xa0", "twopi") + " ; xa0 = psi - a0 ; "
        expr += self._wrapExpression("dm0", "xm0", "twopi") + " ; xm0 = psi - mid0 ; mid0 = (a0 + b0)/2 ; "
        expr += self._dihedralExpression("phi", "xr2", "yr2", "zr2", "xr1", "yr1", "zr1", "0", "0", "0", "xl1", "yl1", "zl1") + " ; "
        expr += self._diffvExpression("xr2", "yr2", "zr2", "x3", "y3", "z3",  "x1", "y1", "z1") + " ; "
        expr += self._diffvExpression("xr1", "yr1", "zr1", "x2", "y2", "z2",  "x1", "y1", "z1") + " ; "
        expr += self._diffvExpression("xl1", "yl1", "zl1", "x5", "y5", "z5",  "x4", "y4", "z4") + " ; "
        expr += "twopi = 2*pi ;"
        expr += "pi = %f" % math.pi
        if kphi is not None:
            if self.CMAnglePhiForce is None:
               self.CMAnglePhiForce = mm.CustomCentroidBondForce(5, expr)
               phiforce = self.CMAnglePhiForce
               numphigroups = 0
               phiforce.addPerBondParameter("kf")
               phiforce.addPerBondParameter("a0")
               phiforce.addPerBondParameter("b0")
               phiforce.setForceGroup(1)
               self.system.addForce(phiforce)
            else:
                phiforce = self.CMAnglePhiForce
                numphigroups = phiforce.getNumGroups()
            phiforce.addGroup(rcpt_cm_groups[0])
            phiforce.addGroup(rcpt_cm_groups[1])
            phiforce.addGroup(rcpt_cm_groups[2])
            phiforce.addGroup( lig_cm_groups[0])
            phiforce.addGroup( lig_cm_groups[1])
            kf = kphi/(kilojoule_per_mole/radians**2)
            a0 = (phi0-phitol)/radians
            b0 = (phi0+phitol)/radians
            groupids = [ numphigroups + i for i in range(5) ]
            phiforce.addBond(groupids, np.array([kf, a0, b0], dtype=np.double) )

        #psi restraint
        # (kf/2) (psi - psi0)^2 with a tolerance and periodicity
        psiforce = None
        expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
        expr += self._wrapExpression("db0", "xb0", "twopi") + " ; xb0 = psi - b0 ; "
        expr += self._wrapExpression("da0", "xa0", "twopi") + " ; xa0 = psi - a0 ; "
        expr += self._wrapExpression("dm0", "xm0", "twopi") + " ; xm0 = psi - mid0 ; mid0 = (a0 + b0)/2 ; "
        expr += self._dihedralExpression("psi", "xr1", "yr1", "zr1", "0", "0", "0", "xl1", "yl1", "zl1", "xl2", "yl2", "zl2") + " ; "
        expr += self._diffvExpression("xr1", "yr1", "zr1", "x2", "y2", "z2",  "x1", "y1", "z1") + " ; "
        expr += self._diffvExpression("xl1", "yl1", "zl1", "x4", "y4", "z4",  "x3", "y3", "z3") + " ; "
        expr += self._diffvExpression("xl2", "yl2", "zl2", "x5", "y5", "z5",  "x3", "y3", "z3") + " ; "
        expr += "twopi = 2*pi ;"
        expr += "pi = %f" % math.pi
        if kpsi is not None:
            if self.CMAnglePsiForce is None:
               self.CMAnglePsiForce = mm.CustomCentroidBondForce(5, expr)
               psiforce = self.CMAnglePsiForce
               numpsigroups = 0
               psiforce.addPerBondParameter("kf")
               psiforce.addPerBondParameter("a0")
               psiforce.addPerBondParameter("b0")
               psiforce.setForceGroup(1)
               self.system.addForce(psiforce)
            else:
                psiforce = self.CMAnglePsiForce
                numpsigroups = psiforce.getNumGroups()
            psiforce.addGroup(rcpt_cm_groups[0])
            psiforce.addGroup(rcpt_cm_groups[1])
            psiforce.addGroup( lig_cm_groups[0])
            psiforce.addGroup( lig_cm_groups[1])
            psiforce.addGroup( lig_cm_groups[2])
            kf = kpsi/(kilojoule_per_mole/radians**2)
            a0 = (psi0-psitol)/radians
            b0 = (psi0+psitol)/radians
            groupids = [ numpsigroups + i for i in range(5) ]
            psiforce.addBond(groupids, np.array([kf, a0, b0], dtype=np.double) )
        return (thetaforce, phiforce, psiforce)

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
            x0 = refpos[p][0]/nanometer
            y0 = refpos[p][1]/nanometer
            z0 = refpos[p][2]/nanometer
            fc1 = fc/(kilojoule_per_mole/nanometer**2)
            tol1 = tol/nanometer
            posrestforce.addParticle(p, np.array([x0, y0, z0, fc1, tol1], dtype=np.double)  )
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
