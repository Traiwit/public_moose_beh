"""some vector operations, vectors are three element lists.

usage:

>>> from bae.vecmath_01 import *
>>> a = [1, 5, 0]; b = [0, 0.5, 0]
>>> # c = a+3*b
>>> c = vector_plus(a, vector_scale(b,3))
"""

__version__="1.18"

_version_history_ = """\
Versions:

1.1 : vector operations as in many of Stephans python scripts
1.2 modif - GP: pow(x,2) replaced by x**2
1.3 added - GP: dist() from findAtMesh
1.4 added - GP: vector_minus
1.5 GP added: in-place operation vector_modif
1.6 GP added: vector_rot_z()
1.7 GP added some matrix operations
1.8 GP added: mat_inverse_2x2, mat_multvec_2x2, mat_transpose
1.9 GP added: vector_sum(), vector_rot_x(), vector_rot_y(), mat_multmat(),
              trans_vert2dir()
       fixed: mat_inverse, vector_rot_y()
1.10 GP added: vector_modif_add_scaled()
1.11 GP added: centre()
1.12 GP added: mat_fromSymTensor()
1.13 GP added trans_rodrigues()
1.14 GP added mat_toRodrigues, mat_orthoNormBaseGS
1.15 FR added mat_multscalar(), mat_submat(), mat_unit()
1.16 FR added functions to calculate Eigenvalues and Eigenvectors
1.17 GP added mat_fromAbqSymTensor, mat_toAbqSymTensor,
        fixed mat_unit, renamed mat_eigenwerte to mat_eigenvalues
1.18 GP added matAbqSym_multvec
"""

## ---------------------------------------------------------------------------
## vector functions v=[x,y,z]

from math import sqrt, sin, cos, pi, acos
from itertools import izip
import math

def dist(A,B):
    l=sqrt((B[0]-A[0])**2+(B[1]-A[1])**2+(B[2]-A[2])**2)
    return l

def vector(A,B):
    """
    @Note: vector(a,b) == vector_minus(b,a)
    """
    return [B[0]-A[0],B[1]-A[1],B[2]-A[2]]

def vector_plus(a,b):
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]]

def vector_minus(a,b):
    """
    @Note: vector_minus(a,b) == vector(b,a)
    """
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]]

def vector_sum(*vecs):
    return [sum(xx) for xx in zip(*vecs)]

def vector_scale(a,r):
    return [a[0]*r,a[1]*r,a[2]*r]

def centre(*vecs):
    """Return the average of given vectors.

    Beware of integer components!
     >>> print centre([0,0,0], [2,2,0], [4,4,2])
     [2,2,0]
     >>> print centre([0,0,0.0], [2,2,0], [4,4,2])
     [2, 2, 0.6666666666666666]
    """
    return [sum(xx)/len(vecs) for xx in zip(*vecs)]

def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def cross(a,b):
    return [-a[2]*b[1]+a[1]*b[2], a[2]*b[0]-a[0]*b[2], -a[1]*b[0]+a[0]*b[1]]

def norm(a):
    l=sqrt(abs(a[0]**2)+abs(a[1]**2)+abs(a[2]**2))
    return [a[0]/l,a[1]/l,a[2]/l]

def length(a):
    l=sqrt(abs(a[0]**2)+abs(a[1]**2)+abs(a[2]**2))
    return l

# in-place vector operation
def vector_modif_add(a, b):
    a[0]+=b[0]; a[1]+=b[1]; a[2]+=b[2]

def vector_modif_scale(a, r):
    a[0]*=r; a[1]*=r; a[2]*=r

def vector_modif_add_scaled(a, b, r):
    "adds r*b to a:  a[i] = a[i]+r*b[i]"
    a[0]+=r*b[0]; a[1]+=r*b[1]; a[2]+=r*b[2]

# rotation
def vector_rot_x(v, deg):
    """rotate vector around the x-axis

    @param v: vector as 3-tuple
    @param deg: rotation angle in degrees
    @returns: turned vector
    """
    alpha = deg*pi/180; c=cos(alpha); s=sin(alpha)
    # rotmat = [[1, 0,  0],
    #           [0, c, -s],
    #           [0, s,  c]]
    return [v[0],
            c*v[1]-s*v[2],
            s*v[1]+c*v[2]]

def vector_rot_y(v, deg):
    """rotate vector around the y-axis

    @param v: vector as 3-tuple
    @param deg: rotation angle in degrees
    @returns: turned vector
    """
    alpha = deg*pi/180; c=cos(alpha); s=sin(alpha)
    # rotmat = [[ c, 0, s],
    #           [ 0, 1, 0],
    #           [-s, 0, c]]
    return [c*v[0]+s*v[2],
            v[1],
            -s*v[0]+c*v[2]]

def vector_rot_z(v, deg):
    """rotate vector around the z-axis

    @param v: vector as 3-tuple
    @param deg: rotation angle in degrees
    @returns: turned vector
    """
    alpha = deg*pi/180; c=cos(alpha); s=sin(alpha)
    # rotmat = [[c, -s, 0],
    #           [s,  c, 0],
    #           [0,  0, 1]]
    return [c*v[0]-s*v[1],
            s*v[0]+c*v[1],
            v[2]]

## ---------------------------------------------------------------------------
## some matrix operations
def mat_submat(mat1,mat2):
    """subtracts a 3x3 from a 3x3 matrix

    @param mat1: a_ij
    @param mat2: b_ij
    @returns: m_ij = a_ij-b_ij
    """
    res=mat_fromSymTensor([0,0,0,0,0,0])
    for i in [0,1,2]:
        for j in [0,1,2]:
            res[i][j]=mat1[i][j]-mat2[i][j]
    return res

def mat_multscalar(mat,scalar):
    """multiplies a 3x3 matrix with a scalar

    @param mat: a_ij
    @param scalar: ...guess...
    @returns: m_ij = scalar*a_ij
    """
    res=mat_fromSymTensor([0,0,0,0,0,0])
    for i in [0,1,2]:
        for j in [0,1,2]:
            res[i][j]=scalar*mat[i][j]
    return res

def mat_transpose(mat):
    """calculate the transpose of a NxM matrix

    @param mat: [[a11,a12,a13,...], [a21,a22,a23,...], ...]
    @returns: the transpose [[a11,a21,...], [a12,a22,...], [a13,a23,...], ...]
    """
    # return zip(*mat)
    return [list(x) for x in izip(*mat)]

def mat_inverse(mat, tol=1E-6):
    """calculate the inverse of a 3x3 matrix

    @param mat: [[a11,a12,a13], [a21,a22,a23], [a31,a32,a33]]
    @param tol: If the absolute value of the determinant of mat is smaller than
    tol, raise ValueError.
    """
    detMat = (mat[0][0]*mat[1][1]*mat[2][2]
            + mat[0][1]*mat[1][2]*mat[2][0]
            + mat[0][2]*mat[1][0]*mat[2][1]
            - mat[2][0]*mat[1][1]*mat[0][2]
            - mat[2][1]*mat[1][2]*mat[0][0]
            - mat[2][2]*mat[1][0]*mat[0][1])
    if abs(detMat) < tol:
        # negative or zero volume, ignore this element
        raise ValueError("Cannot inverse the matrix, the determinant is"
                         " too close to zero: %g" % detMat)

    invMat = [[float(mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1])/detMat,
               float(mat[2][1]*mat[0][2] - mat[2][2]*mat[0][1])/detMat,
               float(mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1])/detMat],
              [float(mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2])/detMat,
               float(mat[2][2]*mat[0][0] - mat[2][0]*mat[0][2])/detMat,
               float(mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2])/detMat],
              [float(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0])/detMat,
               float(mat[2][0]*mat[0][1] - mat[2][1]*mat[0][0])/detMat,
               float(mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0])/detMat]]

    return invMat

def mat_inverse_2x2(mat, tol=1E-6):
    """calculate the inverse of a 2x2 or matrix

    @param mat: [[a11,a12], [a21,a22]]
    @param tol: If the absolute value of the determinant of mat is smaller than
    tol, raise ValueError.
    """
    detMat = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]
    if abs(detMat) < tol:
        # negative or zero volume, ignore this element
        raise ValueError("Cannot inverse the matrix, the determinant is"
                         " too close to zero: %g" % detMat)

    invMat = [[float(mat[1][1])/detMat, -float(mat[0][1])/detMat],
              [-float(mat[1][0])/detMat, float(mat[0][0])/detMat]]

    return invMat

def mat_determinant(mat):
    """calculate the determinant of a 3x3 matrix

    @param mat: [[a11,a12,a13], [a21,a22,a23], [a31,a32,a33]]
    """
    detMat = (mat[0][0]*mat[1][1]*mat[2][2]
            + mat[0][1]*mat[1][2]*mat[2][0]
            + mat[0][2]*mat[1][0]*mat[2][1]
            - mat[2][0]*mat[1][1]*mat[0][2]
            - mat[2][1]*mat[1][2]*mat[0][0]
            - mat[2][2]*mat[1][0]*mat[0][1])
    return detMat

def mat_multvec(mat,vec):
    """performs the matrix product mat * vec;
    res[i] = sum_j mat[i][j]*vec[j];
    mat must be 3x3, vec must be of dim 3
    """
    res = [mat[0][0]*vec[0]+mat[0][1]*vec[1]+mat[0][2]*vec[2],
           mat[1][0]*vec[0]+mat[1][1]*vec[1]+mat[1][2]*vec[2],
           mat[2][0]*vec[0]+mat[2][1]*vec[1]+mat[2][2]*vec[2]]
    return res

def matAbqSym_multvec(matAbqSym, vec):
    """performs the matrix product matAbqSym * vec;
    matAbqSym is a symmetrical tensor in Abaqus/Explicit 6-component
    Voigt-like notation, vec must be of dim 3

    Taken from Vlad's BasicRoutinesJointedMaterial.c PrVe3TeCo6()
    """
    res = [matAbqSym[3]*vec[1] + matAbqSym[5]*vec[2] + matAbqSym[0]*vec[0],
           matAbqSym[1]*vec[1] + matAbqSym[4]*vec[2] + matAbqSym[3]*vec[0],
           matAbqSym[4]*vec[1] + matAbqSym[2]*vec[2] + matAbqSym[5]*vec[0]]
    return res

def mat_multmat(mat1,mat2):
    """performs the matrix product res = mat1 * mat2;
    res[i][j] = sum_k mat1[i][k]*mat2[k][j];
    both mat1 and mat2 must be 3x3
    """
    res = [
        [mat1[0][0]*mat2[0][0]+mat1[0][1]*mat2[1][0]+mat1[0][2]*mat2[2][0],
         mat1[0][0]*mat2[0][1]+mat1[0][1]*mat2[1][1]+mat1[0][2]*mat2[2][1],
         mat1[0][0]*mat2[0][2]+mat1[0][1]*mat2[1][2]+mat1[0][2]*mat2[2][2]],
        [mat1[1][0]*mat2[0][0]+mat1[1][1]*mat2[1][0]+mat1[1][2]*mat2[2][0],
         mat1[1][0]*mat2[0][1]+mat1[1][1]*mat2[1][1]+mat1[1][2]*mat2[2][1],
         mat1[1][0]*mat2[0][2]+mat1[1][1]*mat2[1][2]+mat1[1][2]*mat2[2][2]],
        [mat1[2][0]*mat2[0][0]+mat1[2][1]*mat2[1][0]+mat1[2][2]*mat2[2][0],
         mat1[2][0]*mat2[0][1]+mat1[2][1]*mat2[1][1]+mat1[2][2]*mat2[2][1],
         mat1[2][0]*mat2[0][2]+mat1[2][1]*mat2[1][2]+mat1[2][2]*mat2[2][2]],
        ]
    return res

def mat_multvec_2x2(mat,vec):
    """performs the 2D matrix product mat * vec

    @Note: ignores any components other than with (both) indexes either 0 or 1
    and allways yields a vector with two components.
    """
    res = [mat[0][0]*vec[0]+mat[0][1]*vec[1],
           mat[1][0]*vec[0]+mat[1][1]*vec[1]]
    return res

def mat_linsolve(mat, vec, tol=1E-100):
    """solves res for a 3x3 matrix mat and a component-3 vector vec:
    mat * res = vec
    """

    detMat = (mat[0][0]*mat[1][1]*mat[2][2]
            + mat[0][1]*mat[1][2]*mat[2][0]
            + mat[0][2]*mat[1][0]*mat[2][1]
            - mat[2][0]*mat[1][1]*mat[0][2]
            - mat[2][1]*mat[1][2]*mat[0][0]
            - mat[2][2]*mat[1][0]*mat[0][1])

    if abs(detMat) < tol:
        # negative or zero volume, ignore this element
        raise ValueError("Cannot inverse the matrix, the determinant is too"
                         " close to zero: %g" % detMat)

    inv = [[mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1],
            mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2],
            mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]],
           [mat[2][1]*mat[0][2] - mat[2][2]*mat[0][1],
            mat[2][2]*mat[0][0] - mat[2][0]*mat[0][2],
            mat[2][0]*mat[0][1] - mat[2][1]*mat[0][0]],
           [mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1],
            mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2],
            mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]]]

    res = [float(inv[0][0]*vec[0]+inv[0][1]*vec[1]+inv[0][2]*vec[2])/detMat,
           float(inv[1][0]*vec[0]+inv[1][1]*vec[1]+inv[1][2]*vec[2])/detMat,
           float(inv[2][0]*vec[0]+inv[2][1]*vec[1]+inv[2][2]*vec[2])/detMat]

    return res

def mat_orthoNormBaseGS(mat):
    """Orthogonalize and normalize a 3x3 vector base. Implements the
    Gram-Schmidt process for the first two base vectors. The thrid is taken
    as the cross product of the first two.

    @param mat: List of initial base vectors. The vectors being rows of the
      matrix.
    @returns: List of orthonormal base vectors. I.e. its a matrix with the
      new base vectors as rows.
    """
    e1 = norm(mat[0])
    e2 = norm(vector_minus(mat[1], vector_scale(e1, dot(mat[1], e1))))
    e3 = cross(e1,e2)
    return [e1, e2, e3]


## ---------------------------------------------------------------------------
#{ matrix generation

def mat_fromSymTensor(v):
    """Creates a 3x3 list of lists representing a matrix from a six-component
    vector.

    See also L{mat_fromAbqSymTensor}.

    B{WARNING:} mat_fromSymTensor is not compatible with Abaqus/Explicit tensor
    notation.

    @param v: list of six tensor components in the following order:
       11, 22, 33, 12, 13, 23
    """
    return [ [v[0], v[3], v[4]],
             [v[3], v[1], v[5]],
             [v[4], v[5], v[2]] ]

def mat_fromAbqSymTensor(v):
    """Creates a 3x3 list of lists representing a matrix from a six-component
    vector in the order Abaqus/explicit uses for example in the vumat.

    See also L{mat_fromSymTensor} and the inverse L{mat_toAbqSymTensor}.

    @param v: list of six tensor components in the following order:
       11, 22, 33, 12, 23, 13
    """
    return [ [v[0], v[3], v[5]],
             [v[3], v[1], v[4]],
             [v[5], v[4], v[2]] ]

def mat_toAbqSymTensor(m):
    """Create a six-component vector in the order Abaqus/explicit uses from
    a 3x3 list of lists representing a matrix.
    This is the inverse operation of L{mat_fromAbqSymTensor}.

    @param m: 3x3 matrix (preferably symmetric)
    @returns: list of six tensor components in the following order:
       11, 22, 33, 12, 23, 13
    """
    return [m[0][0], m[1][1], m[2][2],
            0.5*(m[0][1]+m[1][0]),
            0.5*(m[1][2]+m[2][1]),
            0.5*(m[0][2]+m[2][0])]

def mat_unit():
    """Returns the 3x3 unit matrix
    """
    return [ [1.,0.,0.],
             [0.,1.,0.],
             [0.,0.,1.] ]

#}

## ---------------------------------------------------------------------------
#{ Eigenisation of a symmetrix 3x3 matrix

def mat_eigenvalues(A):
    """returns the Eigenvalues of the sym. 3x3 matrix A

    @param A: 3x3 matrix
    @returns: the Eigenvalues in descending order (eig1, eig2, eig3)
         satisfying eig3 <= eig2 <= eig1
    """

    p1 = A[0][1]**2 + A[0][2]**2 + A[1][2]**2
    if (p1 == 0):
        # A is diagonal.
        eig1 = A[0][0]
        eig2 = A[1][1]
        eig3 = A[2][2]
    else:
        q = (A[0][0]+A[1][1]+A[2][2])/3.
        p2 = (A[0][0] - q)**2 + (A[1][1] - q)**2 + (A[2][2] - q)**2 + 2. * p1
        p = math.sqrt(p2 / 6.)
        # (1 / p) * (A - q * I)
        B = mat_fromSymTensor([0,0,0,0,0,0])
        for i in [0,1,2]:
            for j in [0,1,2]:
                if i==j:
                    B[i][j]=(1./p)*(A[i][j]-q)
                else:
                    B[i][j]=(1./p)*(A[i][j])
        detB = mat_determinant(B)
        r = detB / 2.

        # In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        # but computation error can leave it slightly outside this range.
        if (r <= -1.):
            phi = pi / 3.
        elif (r >= 1.):
            phi = 0
        else:
            phi = math.acos(r)/3.

        # the Eigenvalues satisfy eig3 <= eig2 <= eig1
        eig1 = q + 2 * p * math.cos(phi)
        eig3 = q + 2 * p * math.cos(phi + (2*math.pi/3.))
        eig2 = 3 * q - eig1 - eig3     # since trace(A) = eig1 + eig2 + eig3

        return eig1,eig2,eig3

# DEPRECATED alias for mat_eigenvalues
mat_eigenwerte = mat_eigenvalues

def ComputeOrthogonalComplement(W):
    if (abs(W[0])>abs(W[1])):
        # The component of maximum absolute value i s either W[0] or W[ 2 ] .
        invLength = 1. / math.sqrt(W[0]*W[0] + W[2]*W[2])
        U = [-W[2]*invLength,0,+W[0]*invLength]
    else:
        # The component of maximum absolute value i s either W[1] or W[ 2 ]
        invLength = 1. / math.sqrt(W[1]*W[1] + W[2]*W[2])
        U = [0,W[2]*invLength,-W[1]*invLength]

    V = cross(W,U)

    return U,V

def mat_eigenvector(A,eig):
    """Computes the Eigenvector of the symmetrical matrix A corresponding
    to the Eigenvalue eig.

    @param A: symmetrical matrix A
    @param eig: Eigenvalue corresponding to the Eigenvector to be
        computed.

    @returns: the Eigenvector as list of three floats
    """

    row0 = [A[0][0] - eig,A[0][1],A[0][2]]
    row1 = [A[0][1],A[1][1]-eig,A[1][2]]
    row2 = [A[0][2],A[1][2],A[2][2] - eig]
    r0xr1 = cross(row0,row1)
    r0xr2 = cross(row0,row2)
    r1xr2 = cross(row1,row2)
    d0 = dot(r0xr1,r0xr1)
    d1 = dot(r0xr2,r0xr2)
    d2 = dot(r1xr2,r1xr2)
    dmax = d0
    imax = 0
    if (d1>dmax):
        dmax = d1
        imax = 1
    if (d2>dmax):
        imax = 2
    if (imax == 0):
        ev = vector_scale(r0xr1,math.sqrt(d0))
    elif(imax == 1):
        ev = vector_scale(r0xr2,math.sqrt(d1))
    else:
        ev = vector_scale(r1xr2,math.sqrt(d2))

    return ev

def mat_eigenvector1(A, eigenvector0, eigenvalue1):
    """Computes the second Eigenvector of the symmetrical matrix A
    corresponding to the Eigenvalue eig.

    This method is either more efficient or more accurate than
    L{mat_eigenvector}. Fred -the author- should know.

    @param A: symmetrical matrix A
    @param eigenvector0: other Eigenvector already computed earlier
    @param eigenvalue1: Eigenvalue corresponding to the Eigenvector to be
        computed.

    @returns: the Eigenvector as list of three floats
    """

    U,V = ComputeOrthogonalComplement(eigenvector0)

    AU = mat_multvec(A,U)
    AV = mat_multvec(A,V)

    m00 = dot(U,AU)-eigenvalue1
    m01 = dot(U,AV)
    m11 = dot(V,AV)-eigenvalue1

    absM00 = abs(m00)
    absM01 = abs(m01)
    absM11 = abs(m11)


    if (absM00 >= absM11):
        maxAbsComp = max(absM00,absM01)
        if (maxAbsComp > 0):
            if (absM00 >= absM01):
                m01 = m01 / m00
                m00 = 1. / sqrt(1 + m01*m01)
                m01 = m01 * m00
            else:
                m00 = m00 / m01
                m01 = 1 / sqrt(1 + m00*m00)
                m00 = m00 * m01

            eigenvector1 = vector_minus(vector_scale(U,m01),vector_scale(V,m00))
        else:
            eigenvector1 = U
    else:
        maxAbsComp = max(absM11,absM01)
        if (maxAbsComp>0):
            if (absM11>= absM01):
                m01 = m01 / m11
                m11 = 1. / sqrt(1 + m01*m01)
                m01 = m01 * m11
            else:
                m11 = m11 / m01
                m01 = 1. / sqrt(1 + m11*m11)
                m11 = m11 * m01
            eigenvector1 = vector_minus(vector_scale(U,m11),vector_scale(V,m01))
        else:
            eigenvector1 = U

    return eigenvector1

def mat_eigenisation1(A):
    """returns the Eigenvalues and Eigenvectors of the sym. 3x3 matrix A

    Usage:
     >>> eig1,eig2,eig3,ev1,ev2,ev3 = mat_eigenisation1(A)

    @param A: 3x3 matrix
    @returns: eig1,eig2,eig3,ev1,ev2,ev3
              the Eigenvalues satisfy eig3 <= eig2 <= eig1 and the
              corresponding Eigenvectors ev1, ev2, ev3
    """
    eig1,eig2,eig3 = mat_eigenvalues(A)

    ev1=norm(mat_eigenvector(A,eig1))
    ev2=norm(mat_eigenvector(A,eig2))
    ev3=norm(mat_eigenvector(A,eig3))

    return eig1,eig2,eig3,ev1,ev2,ev3

def mat_eigenisation2(A):
    """Returns the Eigenvalues and Eigenvectors of the sym. 3x3 matrix A

    Usage:
     >>> eig1,eig2,eig3,ev1,ev2,ev3 = mat_eigenisation1(A)

    Contrary to L{mat_eigenisation1} this function uses L{mat_eigenvector1}
    for the second Eigenvector and uses the cross product of the first two as
    the third.

    @param A: 3x3 matrix
    @returns: eig1,eig2,eig3,ev1,ev2,ev3
              the Eigenvalues satisfy eig3 <= eig2 <= eig1 and the
              corresponding Eigenvectors ev1, ev2, ev3
    """

    eig1,eig2,eig3 = mat_eigenvalues(A)

    ev1 = norm(mat_eigenvector(A,eig1))
    ev2 = mat_eigenvector1(A,ev1,eig2)
    ev3 = cross(ev1,ev2)

    return eig1,eig2,eig3,ev1,ev2,ev3

#}

## ---------------------------------------------------------------------------
#{ special transformation matrices

def trans_vert2dir(direction):
    """'TRANSformation from VERTical TO other DIRection'
    Returns the transformation matrix for a rotation from vertically upwards
    to the given direction. I.e. the z axis is rotated to the given direction.

    The rotation axis is perpendicular to the z axis and perpendicular to
    direction. If direction points to the neg. z-direction (i.e. direction[2]
    == -length(direction)) then the rotaion is performed around the y-axis.

    The i-th column (!) vector A[:][i] of the resulting transformation Matrix A
    is the rotated i-th base vector e_i given in the global (not rotated)
    coordinate system.

    A given vector v is rotated by multiplying the vector to this
    transformation vector:
     >>> #--- bolt initialization
     >>> bolt_axis = [0,0,1]     # initial direction: vertical
     >>> bolt_normal = [1,0,0]
     >>> #--- rotation to final position
     >>> bolt_dir = [0.5, 0, 1]  # new direction
     >>> rotmat = trans_vert2dir(bolt_dir)
     >>> bolt_axis = mat_multvec(rotmat, bolt_axis)
     >>> bolt_normal = mat_multvec(rotmat, bolt_normal)
    """
    d=norm(direction)
    L2 = d[0]**2 + d[1]**2

    if L2<1E-12:
        if d[2]>0.999:
            # d == z-axis
            A = [[1,0,0], [0,1,0], [0,0,1]]
        else:
            # d == -z-axis => rotation around y-axis
            A = [[-1,0,0], [0,1,0], [0,0,-1]]
    else:
        A = [ [(d[1]**2 + d[0]**2*d[2])/L2, (d[2]-1.0)*d[0]*d[1]/L2, d[0] ],
              [(d[2]-1.0)*d[0]*d[1]/L2, (d[0]**2 + d[1]**2*d[2])/L2, d[1] ],
              [-d[0], -d[1], d[2] ] ]
    return A

def trans_rodrigues(direction, deg):
    """'TRANSformation matrix of rotation in Rodrigues formulation'
    see: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula

    Returns the transformation matrix for a rotation given in Rodrigues
    formulation: Rotation of angle alpha around unit vector k (right hand rule)

    cross product matrix K (for any vector v: Kv = cross(k, v)::
          [   0    -k_3    k_2 ]
      K = [  k_3     0    -k_1 ]
          [ -k_2    k_1     0  ]

    Transformation matrix of the rotation is given as
    R = I + sin(alpha) K + (1-cos alpha) K^2

    The i-th column (!) vector R[:][i] of the resulting transformation Matrix R
    is the rotated i-th base vector e_i given in the global (not rotated)
    coordinate system.

    A given vector v is rotated by multiplying the vector to this
    transformation vector:
     >>> rotmat = trans_rodrigues([1,0,0], 90)  # 90 deg rotation around x
     >>> x = [0,0,1]     # initial direction: vertical
     >>> x_rot = mat_multvec(rotmat, x)
     >>> print x_rot
     [0,-1,0]

    @Note: For the inverse of this operation see L{mat_toRodrigues}
    """

    k=norm(direction)
    alpha = deg*pi/180; c=(1.0-cos(alpha)); s=sin(alpha)

    #       [ -k_2^2 - k_3^2       k_1 k_2          k_1 k_3    ]
    # K^2 = [     k_1 k_2      -k_1^2 - k_3^2       k_2 k_3    ]
    #       [     k_1 k_3          k_2 k_3      -k_1^2 - k_2^2 ]
    #
    # R = I + sin(alpha) K + (1-cos alpha) K^2

    R = [[1-c*(k[1]**2+k[2]**2), -s*k[2] + c*k[0]*k[1],  s*k[1] + c*k[0]*k[2]],
         [ s*k[2] + c*k[0]*k[1], 1-c*(k[0]**2+k[2]**2), -s*k[0] + c*k[1]*k[2]],
         [-s*k[1] + c*k[0]*k[2],  s*k[0] + c*k[1]*k[2], 1-c*(k[0]**2+k[1]**2)]]
    return R

def mat_toRodrigues(R):
    """This is the inverse operation of L{trans_rodrigues}: return the rotation
    unit vector and rotation angle of a rotation given by its transformation
    matrix R. No check is being made that R actually performs a rotation.

    Given the transformation matrix R for the rotation we get the rotation
    angle alpha and the unit vector of the rotaion axis k like so:

    alpha = arccos( 1/2 (trace(R)-1) )

    ... and::
                           [ R_32 - R_23 ]
      k = 1/(2*sin(alpha)) [ R_13 - R_31 ]
                           [ R_21 - R_12 ]

    @returns: a (direction, deg)-tuple. direction is a unit vector, deg the
    rotation angle in degree.
    """
    alpha = acos( 0.5*(R[0][0] + R[1][1] + R[2][2] - 1) )
    k = [R[2][1]-R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]]
    if any(x!=0 for x in k):
        k = norm(k)
    return (k, 180*alpha/pi)

#}

## ---------------------------------------------------------------------------
#{ matrix output

def formatMats(*args, **kwargs):
    """Format a multiline string with matrices

    @param args: What to print: Matrices, or other strings which go in the
      centre line. Anything but lists and strings will be ignored. No empty
      list!

    @keyword format: format string for number conversion, defaults to "%8.2g"
    """
    try:
        fmt = kwargs["format"]
    except KeyError:
        fmt="%8.2g"

    args2 = list(args)
    for i, val in enumerate(args):
        if isinstance(val, basestring):
            space = " "*len(val)
            args2[i] = [space, val, space]

    result = ""
    for line in izip(*args2):
        for item in line:
            if isinstance(item, (list, tuple)):
                result += "[ %s ]" % (",".join(fmt % x for x in item))
            elif isinstance(item, (float, int)):
                result += "[ %s ]" % (fmt % item)
            else:
                result += "%s" % item
        result += "\n"
    return result
