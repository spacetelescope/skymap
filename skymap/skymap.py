import numpy as np

def y2theta(yp):
    # yp = yi/pi
    try:
        idx = np.abs(yp) <= 1/4
        sintheta = np.zeros(len(yp))
        sintheta[idx] = 8 * yp[idx]/3
        sintheta[~idx] = (1 - (2 - 4 * np.abs(yp[~idx]))**2/3) * np.abs(yp[~idx])/yp[~idx]
    except:
        if np.abs(yp) <= 1/4:
            sintheta = 8*yp/3
        else:
            sintheta = (1 - (2 - 4 * np.abs(yp))**2/3) * np.abs(yp)/yp
    return np.arcsin(sintheta) * 180/np.pi

def count(N):
    cnt = np.arange(1,N,dtype=np.intc) * 4
    cnt = np.concatenate([cnt, np.ones(2*N+1, dtype=np.intc)*4*N,cnt[::-1]], dtype=np.intc)
    return(cnt)

def ang2point(theta, phi, radius):
    # Polar to cartesian coordinates
    xi = radius * np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(phi))
    yi = radius * np.cos(np.deg2rad(theta)) * np.sin(np.deg2rad(phi))
    zi = radius * np.sin(np.deg2rad(theta))
    points = np.vstack((xi,yi,zi)).T
    return points


def romantessellation(NSIDE, radius=1):
    # Definition of Roman tessellation
    
    yup = (2 * NSIDE - np.arange(0, 4*NSIDE+1)) / (4 * NSIDE) # yi/pi
    thetau = y2theta(yup)
    cnt = count(NSIDE)
    yuvp = (2 * NSIDE - np.arange(0, 4*NSIDE) - 0.5) / (4 * NSIDE) # yi/pi # vertices
    thetav = y2theta(yuvp)
   
    vertices = []
    # First tile is the polar cap
    nphi = [225]
    ntheta = [90]
    thu = thetav[0]
    _theta = np.array([thu] * 36)
    _phi = np.arange(0,360,10)
    vertices.append(ang2point(_theta,_phi,radius))
    
    cra, cdec = [0], [90]
    ramin, ramax = [0], [360]
    decmin, decmax = [thu], [90]
    
    npoints = 20
    for cnt_, thu_, thd, yup_ in zip(cnt,thetau[1:],thetav[1:], yup[1:]):
        _theta = thu_
        dphi = 360 / (2 * cnt_)
        for i in range(2*cnt_):
            dphi = 360 / (2 * cnt_)
            _phi = i * dphi
            cra.append(i*dphi)
            cdec.append(thu_)
            ramin.append(-dphi/2+i*dphi)
            ramax.append(+dphi/2+i*dphi)
            decmin.append(thu)
            decmax.append(thd)
            nphi.append(_phi)
            ntheta.append(_theta)
            vphi = np.arange(-npoints//2,npoints//2+1) * dphi/npoints + i*dphi
            vphi = np.concatenate([vphi,vphi[::-1]])
            vtheta = np.concatenate([[thd]*(npoints+1), [thu]*(npoints+1)])
            vtheta = np.array(vtheta)
            vphi = np.array(vphi)
            #print(i,n,vphi,vtheta)
            vertices.append(ang2point(vtheta, vphi, radius))
        thu = thd
        
    # Last tile: opposite polar cup
    nphi.append(135)
    ntheta.append(-90)
    ntheta = 90 - np.array(ntheta)        
    nphi = np.array(nphi)
    ntheta = np.array(ntheta)
    _theta = np.array([thu] * 36)
    _phi = np.arange(0,360,10)
    vertices.append(ang2point(_theta,_phi,radius))
    cra.append(0)
    cdec.append(-90)
    ramin.append(0)
    ramax.append(0)
    decmin.append(-90)
    decmax.append(thu)
    ramin, ramax = np.array(ramin), np.array(ramax)
    decmin, decmax = np.array(decmin), np.array(decmax)
   
    
    # Save data as a FITS table
    from astropy.io import fits
    col1 = fits.Column(name='cra', format='E', array=cra)
    col2 = fits.Column(name='cdec', format='E', array=cdec)
    col3 = fits.Column(name='ramin', format='E', array=ramin)
    col4 = fits.Column(name='ramax', format='E', array=ramax)
    col5 = fits.Column(name='decmin', format='E', array=decmin)
    col6 = fits.Column(name='decmax', format='E', array=decmax)
    coldefs = fits.ColDefs([col1, col2, col3, col4, col5, col6])
    hdu = fits.BinTableHDU.from_columns(coldefs)
    hdu.writeto('romantessellation.fits', overwrite='True')

    return ntheta, nphi, ramin, ramax, decmin, decmax, vertices

def doublepixelization(NSIDE, radius=1):
    # Double pixelization vertices for Healpix H=4

    yup = (2 * NSIDE - np.arange(0, 4*NSIDE+1)) / (4 * NSIDE) # yi/pi
    thetau = y2theta(yup)
    cnt = count(NSIDE)
    yuvp = (2 * NSIDE - np.arange(0, 4*NSIDE) - 0.5) / (4 * NSIDE) # yi/pi # vertices
    thetav = y2theta(yuvp)

    vertices = []
    nphi = [225]
    ntheta = [90]
    thu = thetav[0]
    _theta = np.array([thu] * 4)
    _phi = np.array([0,90,180,270])
    vertices.append(ang2point(_theta,_phi,radius))
    # transition angle ....
    # thetaX = np.arcsin(2/3)*180/np.pi

    for cnt_, thu_, thd, yup_ in zip(cnt,thetau[1:],thetav[1:], yup[1:]):
        _theta = thu_
        dphi = 360 / (2 * cnt_)
        iup = 0
        idown = 0
        for i in range(2*cnt_):
            dphi = 360 / (2 * cnt_)
            _phi = i * dphi
            nphi.append(_phi)
            ntheta.append(_theta)
            n = cnt_/2 # number of a side length (longest)
            if np.abs(np.abs(yup_) - 0.25) < 0.0001:
                if _theta > 0:
                    dphi = 360 / (2 * cnt_) # 4 is the triangles
                    dphiu = 360 / (2 * cnt_ - 4)
                    if i%n == 0:
                        vphi = [-dphi/2+_phi,_phi,dphi/2+_phi,_phi]
                        vtheta = [thd, thd, thd, thu]
                        idown += 1
                    else:
                        vphi = [-dphi/2+idown*dphi, dphi/2+idown*dphi, (iup+1)*dphiu, iup*dphiu]
                        vtheta = [thd, thd, thu, thu]
                        idown += 1
                        iup += 1
                else:
                    dphi = 360 / (2 * cnt_) # 4 is the triangles
                    dphiu = 360 / (2 * cnt_ - 4)
                    if i%n == 0:
                        vphi = [-dphi/2+_phi,_phi,dphi/2+_phi,_phi]
                        vtheta = [thu,thu,thu,thd]
                        idown += 1
                    else:
                        vphi = [-dphi/2+idown*dphi, dphi/2+idown*dphi, (iup+1)*dphiu, iup*dphiu]
                        vtheta = [thu, thu, thd, thd]
                        idown += 1
                        iup += 1
            elif  yup_ > 0.25:
                dphi = 360 / (2 * cnt_ + 4) # 4 is the triangles
                dphiu = 360 / (2 * cnt_ - 4)
                if i%n == 0:
                    vphi = [-dphi+_phi,_phi,dphi+_phi,_phi]
                    vtheta = [thd, thd, thd, thu]
                    if i == 0:
                        idown += 1
                    else:
                        idown += 2
                else:
                    vphi = [dphi+(idown-1)*dphi, dphi+idown*dphi, (iup+1)*dphiu, iup*dphiu]
                    vtheta = [thd, thd, thu, thu]
                    idown += 1
                    iup += 1
            elif yup_ < -0.25:
                dphi = 360 / (2 * cnt_ + 4) # 4 is the triangles
                dphiu = 360 / (2 * cnt_ - 4)
                if i%n == 0:
                    vphi = [-dphi+_phi,_phi,dphi+_phi,_phi]
                    vtheta = [thu, thu, thu, thd]
                    if i == 0:
                        idown += 1
                    else:
                        idown += 2
                else:
                    vphi = [dphi+(idown-1)*dphi, dphi+idown*dphi, (iup+1)*dphiu, iup*dphiu]
                    vtheta = [thu, thu, thd, thd]
                    idown +=1
                    iup += 1
            else:
                dphi = 360 / (2 * cnt_) 
                vphi = [-dphi/2 + _phi, -dphi/2+_phi+dphi, -dphi/2+_phi+dphi, -dphi/2 + _phi]
                vtheta = [thu, thu, thd, thd]
            vtheta = np.array(vtheta)
            vphi = np.array(vphi)
            #print(i,n,vphi,vtheta)

            #id = vphi < 0
            #vphi[id] = 360 - vphi[id]
            vertices.append(ang2point(vtheta, vphi, radius))
        thu = thd
        
    nphi.append(135)
    ntheta.append(-90)
    ntheta = 90 - np.array(ntheta)        
    nphi = np.array(nphi)
    ntheta = np.array(ntheta)
    _theta = np.array([thu] * 4)
    _phi = np.array([0,90,180,270])
    vertices.append(ang2point(_theta,_phi,radius))

    return ntheta, nphi, vertices

# Polyhedra

def octahedron(r):
    import numpy as np
    verts = [
        [1,0,0],
        [0,1,0],
        [0,0,1],
        [0,-1,0],
        [0,0,-1],
        [-1,0,0]
    ]
    verts = np.array(verts)
    faces = [ 
        [(1,0,0), (0,1,0), (0,0,1)], 
        [(1,0,0), (0,0,-1), (0,1,0)], 
        [(1,0,0), (0,0,1), (0,-1,0)], 
        [(1,0,0), (0,-1,0), (0,0,-1)], 
        [(-1,0,0), (0,0,1), (0,1,0)], 
        [(-1,0,0), (0,1,0), (0,0,-1)], 
        [(-1,0,0), (0,-1,0), (0,0,1)], 
        [(-1,0,0), (0,0,-1), (0,-1,0)], 
    ]
    
    #faces = []
    #for i in [5,0]:
    #    faces.append([verts[i],verts[2],verts[3]])
    #    faces.append([verts[i],verts[3],verts[4]])
    #    faces.append([verts[i],verts[4],verts[1]])
    #    faces.append([verts[i],verts[1],verts[2]])
    faces = np.array(faces)
    return verts, faces

def octasphere(recursion):
    import numpy as np
    verts, triangles = octahedron(1)
    
    def middlepoint(p1, p2):
        mp = np.zeros(3)
        for i in range(3):
            mp[i] = (p1[i]+p2[i])*0.5
    
        nmp = np.sqrt(np.sum(mp**2))
        mp = mp/nmp
        return mp
    
    for i in range(recursion):
        faces = []
        for triangle in triangles:
            a = middlepoint(triangle[0], triangle[1])
            b = middlepoint(triangle[1], triangle[2])
            c = middlepoint(triangle[2], triangle[0])
              
            faces.append([triangle[0], a, c])
            faces.append([triangle[1], b, a])
            faces.append([triangle[2], c, b])
            faces.append([a, b, c])
        triangles = faces
        
    return triangles


# http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html

def icosahedron(r):
    import numpy as np
    # (  0, ±1, ±φ )
    # ( ±1, ±φ,  0 )
    # ( ±φ,  0, ±1 )
    # φ = (1 + √5) / 2
    phi = (1+np.sqrt(5))/2 
    verts = []
    for i in [1,-1]:
        for j in [-1,1]:
            verts.append(r * [j, i* phi,0])
    for i in [1,-1]:
        for j in [-1,1]:
            verts.append(r * [0,j,i*phi])
    for i in [1,-1]:
        for j in [-1,1]:            
            verts.append(r * [i* phi, 0, j])
    verts = np.array(verts)
    norm = np.sqrt(1 + phi**2)
    verts /= norm
    
    faces = []
    faces.append([verts[0], verts[11], verts[5]])
    faces.append([verts[0], verts[5], verts[1]])
    faces.append([verts[0], verts[1], verts[7]])
    faces.append([verts[0], verts[7], verts[10]])
    faces.append([verts[0], verts[10], verts[11]])

    faces.append([verts[1], verts[5], verts[9]])
    faces.append([verts[5], verts[11], verts[4]])
    faces.append([verts[11], verts[10], verts[2]])
    faces.append([verts[10], verts[7], verts[6]])
    faces.append([verts[7], verts[1], verts[8]])

    faces.append([verts[3], verts[9], verts[4]])
    faces.append([verts[3], verts[4], verts[2]])
    faces.append([verts[3], verts[2], verts[6]])
    faces.append([verts[3], verts[6], verts[8]])
    faces.append([verts[3], verts[8], verts[9]])
                  
    faces.append([verts[4], verts[9], verts[5]])
    faces.append([verts[2], verts[4], verts[11]])
    faces.append([verts[6], verts[2], verts[10]])
    faces.append([verts[8], verts[6], verts[7]])
    faces.append([verts[9], verts[8], verts[1]])

    faces = np.array(faces)
    return verts, faces


# Icosphere

def icosphere(recursion):
    import numpy as np
    verts, triangles = icosahedron(1)
    
    def middlepoint(p1, p2):
        mp = np.zeros(3)
        for i in range(3):
            mp[i] = (p1[i]+p2[i])*0.5
    
        nmp = np.sqrt(np.sum(mp**2))
        mp = mp/nmp
        return mp
    
    for i in range(recursion):
        faces = []
        for triangle in triangles:
            a = middlepoint(triangle[0], triangle[1])
            b = middlepoint(triangle[1], triangle[2])
            c = middlepoint(triangle[2], triangle[0])
              
            faces.append([triangle[0], a, c])
            faces.append([triangle[1], b, a])
            faces.append([triangle[2], c, b])
            faces.append([a, b, c])
        triangles = faces
        
    return triangles
