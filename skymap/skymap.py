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
    ramax.append(360)
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


def tiles2asdf(theta, phi, ramin, ramax, decmin, decmax,pixsize=0.055, cellsize=4800, border=100, outfile='skymap.asdf'):
    """
    Generate an asdf file with the metadata of all the skycell files for all the tiles.
    Input:
         theta        latitude of tile
         phi          longitude of tile
         ramin, ramax  limits in RA of tile
         decmin, decmax  limits in Declination of tile
         pixsize       pixel size in arcsec (default is 0.055 arcsec)
         cellsize      sky cell size in pixels (default is 4800)
         border        border of sky cell overlapping adjacent cells in pixels (default is 100)
         outfile       name of ASDF output file 
    """
    from astropy.wcs import WCS
    import asdf
    import numpy as np
    from datetime import time,date,datetime
    
    d = date(2025, 1, 1)
    t = time(0, 0, 0)
    startdate = datetime.combine(d, t)

    # Structure of the tile metadata
    wcs_tile_dtype = [
        ('index', 'i4'),
        ('ra_tangent', 'f8'),      # RA tile center  [crval1]
        ('dec_tangent', 'f8'),     # Dec tile center [crval2]
        ('ra_min', 'f8'),         # Limits in RA of the tile
        ('ra_max', 'f8'),         # 
        ('dec_min', 'f8'),        # Limits in Dec of the tile
        ('dec_max', 'f8'),        #
        ('orientat', 'f4'),       # Orientation projection [crota2]
        ('x_tangent', 'f8'), # x coord of projection region center
        ('y_tangent', 'f8'), # y coord of projection region center
        ('nx', 'i4'), # x-size of projection region in pixels
        ('ny', 'i4'), # y-size of projection region in pixels
        ('skycell_start','i4'),           # First cell index
        ('skycell_end','i4'),             # Last cell index
        #('nxy_skycell','i4'),            # square cell side in pixels [naxis1 and naxis2]
        #('skycell_border_pixels', 'i4'),         # border of cell in pixels
        ('pixel_scale','f4')           # pixel scale in degrees
    ]

    # Structure of the cell metadata
    wcs_cell_dtype = [
        ('name', 'U16'),
        ('ra_center', 'f8'),     # center of the cell
        ('dec_center', 'f8'),    # center of the cell
        ('orientat', 'f4'),      # orientation of the cell
        ('x_tangent', 'f8'), # x coordinate of projection center in pixels [crpix1]
        ('y_tangent', 'f8'), # y coordinate of projection center in pixels [crpix2]
        ('ra_corn1', 'f8'),      # RA corner 1
        ('dec_corn1', 'f8'),
        ('ra_corn2', 'f8'),
        ('dec_corn2', 'f8'),
        ('ra_corn3', 'f8'),
        ('dec_corn3', 'f8'),
        ('ra_corn4', 'f8'),
        ('dec_corn4', 'f8')
    ]

    wcs = WCS(naxis=2)
    pix = pixsize/3600 # Pixel size
    wcs.wcs.cdelt = [pix,pix]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    ntiles = len(theta)
    hcellsize = cellsize / 2
    nxy = cellsize + border * 2
    # Grid of cells
    n = 35  # enough to cover typical N=13 tile size with nx=4800
    row, col = np.indices((2*n+1, 2*n+1))
    row, col = row - n, col - n
    # Name format
    #namefmt = 'a{0:03d}d{1:s}{2:02d}x{3:s}{4:02d}y{5:s}{6:02d}'
    namefmt = '{0:03d}{1:s}{2:02d}x{3:02d}y{4:02d}'
    # Tiles and cells structured numpy arrays
    tiles = np.empty(ntiles, dtype=wcs_tile_dtype)
    cells = np.empty(0,  dtype=wcs_cell_dtype)
    for itile in range(ntiles):
        tile = tiles[itile]
        # Counter to check the progress ...
        if (itile % 1000) == 0:
            print('\n'+str(itile), end='')
        elif (itile %100) == 0:
            print(':', end='')
        elif (itile %25) == 0:
            print('.', end='')  
        ra0, dec0 = phi[itile], 90-theta[itile]
        ramin_, ramax_ = ramin[itile], ramax[itile]
        decmin_, decmax_ = decmin[itile], decmax[itile]
        if decmin_ > decmax_: decmin_, decmax_ = decmax_, decmin_
        # Central pixel
        dec_ = min(np.abs(decmin_), np.abs(decmax_))
        if itile in [0, ntiles-1]:
            ny = np.abs(decmin_ - decmax_) / pix * 2
            nx = ny
        else:
            ny = np.abs(decmin_ - decmax_) / pix
            nx = np.abs(ramax_ - ramin_) * np.cos(dec_ * np.pi/180) / pix
        # Make nx and ny even number
        nx = nx // 2 * 2
        ny = ny // 2 * 2
        x0t, y0t = (nx - 1) / 2, (ny - 1) / 2
        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = [pix,pix]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crval = [ra0, dec0]
        wcs.wcs.crpix = [x0t, y0t] # Between two central pixels
        wcs.array_shape = [nx,ny]

        #n = 35
        #row, col = np.indices((2*n+1, 2*n+1))
        #row, col = row - n, col - n
        x0, y0 = row * cellsize + x0t, col * cellsize + y0t

        # Cell corners
        x1, y1 = x0 - hcellsize, y0 - hcellsize
        x2, y2 = x0 + hcellsize, y0 - hcellsize
        x3, y3 = x0 + hcellsize, y0 + hcellsize
        x4, y4 = x0 - hcellsize, y0 + hcellsize
        # Extended cell corners (including the overlap with contiguous cells)
        x1e, y1e = x1 - border, y1 - border
        x2e, y2e = x2 + border, y2 - border
        x3e, y3e = x3 + border, y3 + border
        x4e, y4e = x4 - border, y4 + border
        # Cell corners
        a0, d0 = wcs.wcs_pix2world(x0, y0, 0)
        a1, d1 = wcs.wcs_pix2world(x1, y1, 0)
        a2, d2 = wcs.wcs_pix2world(x2, y2, 0)
        a3, d3 = wcs.wcs_pix2world(x3, y3, 0)
        a4, d4 = wcs.wcs_pix2world(x4, y4, 0)
        a1e, d1e = wcs.wcs_pix2world(x1e, y1e, 0)
        a2e, d2e = wcs.wcs_pix2world(x2e, y2e, 0)
        a3e, d3e = wcs.wcs_pix2world(x3e, y3e, 0)
        a4e, d4e = wcs.wcs_pix2world(x4e, y4e, 0)

        if ramin_ < 0:
            idn = a0 > 270
            a0[idn] -= 360
            a1[idn] -= 360
            a2[idn] -= 360
            a3[idn] -= 360
            a4[idn] -= 360
        # Case of pole tiles
        if itile == 0:
            c1 = d1 >= decmin_
            c2 = d2 >= decmin_
            c3 = d3 >= decmin_
            c4 = d4 >= decmin_
        elif itile == ntiles-1:
            c1 = d1 < decmax_
            c2 = d2 < decmax_
            c3 = d3 < decmax_
            c4 = d4 < decmax_
        # Case of generic ring tiles
        else:
            c1 = (a1 >= ramin_) & (a1 < ramax_) & (d1 >= decmin_) & (d1 < decmax_)
            c2 = (a2 >= ramin_) & (a2 < ramax_) & (d2 >= decmin_) & (d2 < decmax_)
            c3 = (a3 >= ramin_) & (a3 < ramax_) & (d3 >= decmin_) & (d3 < decmax_)
            c4 = (a4 >= ramin_) & (a4 < ramax_) & (d4 >= decmin_) & (d4 < decmax_)

        idx, idy = np.where(c1 | c2 | c3 | c4)
        # Part of the file name (tile center coords)
        if dec0 < 0:
            dsign = 'm'
        else:
            dsign = 'p'
        ra0_, dec0_ = round(ra0), round(np.abs(dec0))
        # Fields for tiles
        tile['index'] = itile
        tile['ra_tangent'] = '{0:.10f}'.format(ra0)
        tile['dec_tangent'] = '{0:.10f}'.format(dec0)
        tile['ra_min'] = '{0:.10f}'.format(ramin_)
        tile['ra_max'] = '{0:.10f}'.format(ramax_)
        tile['dec_min'] = '{0:.10f}'.format(decmin_)
        tile['dec_max'] = '{0:.10f}'.format(decmax_)       
        tile['skycell_start'] = len(cells)     # First cell index
        #tile['nxy_skycell'] = nxy
        #tile['skycell_border_pixels'] = border
        tile['pixel_scale'] = pix
        tile['x_tangent'] = x0t
        tile['y_tangent'] = y0t
        tile['nx'] = nx
        tile['ny'] = ny     
        tile['orientat'] = 0 # Orientation projection
        # Generate a numpy structured array for the cells of the tile 
        cell = np.empty(len(idx), dtype=wcs_cell_dtype)
        for icell, (idx_, idy_) in enumerate(zip(idx, idy)):
            #xsign, ysign = 'p', 'p'
            #if x0[idx_, idy_] < 0:
            #    xsign = 'm'
            #if y0[idx_, idy_] < 0:
            #    ysign = 'm'            
            #x0_, y0_ = np.abs(row[idx_, idy_]), np.abs(col[idx_, idy_])
            # Let's assume that (50,50) is the coordinates of the central cell
            # This is done to avoid signs in the x,y coordinates of a cell
            x0_, y0_ = 50+row[idx_, idy_], 50+col[idx_, idy_]
            #cell[icell]['name'] = namefmt.format(ra0_, dsign, dec0_, xsign, x0_, ysign, y0_)
            cell[icell]['name'] = namefmt.format(ra0_, dsign, dec0_, x0_, y0_)
            cell[icell]['ra_center'] = '{0:.10f}'.format(a0[idx_, idy_]) # cell center
            cell[icell]['dec_center'] = '{0:.10f}'.format(d0[idx_, idy_])
            cell[icell]['orientat'] = ra0 - a0[idx_, idy_] # orientation wrt tile
            xpix = x0t - x0[idx_, idy_] + 2499.5
            ypix = y0t - y0[idx_, idy_] + 2499.5
            cell[icell]['x_tangent'] = xpix # position of tile center
            cell[icell]['y_tangent'] = ypix
            cell[icell]['ra_corn1'] = '{0:.10f}'.format(a1e[idx_, idy_])
            cell[icell]['dec_corn1'] = '{0:.10f}'.format(d1e[idx_, idy_])
            cell[icell]['ra_corn2'] = '{0:.10f}'.format(a2e[idx_, idy_])
            cell[icell]['dec_corn2'] = '{0:.10f}'.format(d2e[idx_, idy_])
            cell[icell]['ra_corn3'] = '{0:.10f}'.format(a3e[idx_, idy_])
            cell[icell]['dec_corn3'] = '{0:.10f}'.format(d3e[idx_, idy_])
            cell[icell]['ra_corn4'] = '{0:.10f}'.format(a4e[idx_, idy_])
            cell[icell]['dec_corn4'] = '{0:.10f}'.format(d4e[idx_, idy_])
            # Concatenate cells to the cells from previous tiles
        cells = np.concatenate([cells, cell])
        tile['skycell_end'] = len(cells)       # Last cell index

    # Save the file
    #tree = {}
    roman = asdf.tagged.TaggedDict()
    instrument = {'name':'WFI'}
    meta = {'author': "Dario Fadda",
            'description':"Skycells covering the celestial sphere",
            'homepage': "http://github.com/spacetelescope/skymap",
            'instrument': instrument,
            'nxy_skycell': nxy,
            'origin': 'STSCI',
            'pedigree':'GROUND',
            'plate_scale': 0.55,
            'reftype': 'SKYCELLS',
            'skycell_border_pixels': border,
            'telescope': 'ROMAN',
            'useafter': startdate,
            }
    roman['meta'] = meta
    roman['projection_regions'] = tiles
    roman['skycells'] = cells
    roman['datamodel_name'] = "RomanSkycellsRefModel"
    tree = {
    "roman": roman,
    }
    #tree.update({'projection_regions': tiles})
    #tree.update({'skycells': cells})
    ff = asdf.AsdfFile(tree)
    ff.write_to(outfile)
    return 1

def plotskytile(skymap, itile, distortion=False):
    """
    Plot a skytile (given the index) with all the skycells.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    skytiles = skymap['projection_regions']
    skycells = skymap['skycells']
    tile = skytiles[itile]
    cells = skycells[tile['skycell_start']:tile['skycell_end']]
    c0 = tile['ra_tangent']

    if itile in [0,len(skytiles)-1]:
        dd = 0.04
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        for c in cells:
            cosfac = np.cos(c['dec_corn1'] * np.pi/180)
            c1 = c['ra_corn1']
            c2 = c['ra_corn2']
            c3 = c['ra_corn3']
            c4 = c['ra_corn4']
            if itile == 0:
                d1 = 90-c['dec_corn1']
                d2 = 90-c['dec_corn2']
                d3 = 90-c['dec_corn3']
                d4 = 90-c['dec_corn4']
            else:
                d1 = c['dec_corn1']
                d2 = c['dec_corn2']
                d3 = c['dec_corn3']
                d4 = c['dec_corn4']
            ax.plot(np.array([c1, c2, c3, c4, c1])* np.pi/180,
                    [d1, d2, d3, d4, d1], color='skyblue')
            ac = c['ra_center']* np.pi/180
            dc = 90-c['dec_center']
            orient = (90 - c['orientat']) * np.pi/180
        angle = np.arange(0,2*np.pi,0.01)
        print('Dec max,min is: ', tile['dec_max'], tile['dec_min'])
        if itile == 0:
            dec = 90-tile['dec_min']
        else:
            dec = tile['dec_max']
        ax.plot(angle, dec*np.ones(len(angle)),color='red')
    else:
        fig, ax = plt.subplots()
        for c in cells:
            cosfac = np.cos(c['dec_corn1'] * np.pi/180)
            c1 = c['ra_corn1']
            c2 = c['ra_corn2']
            c3 = c['ra_corn3']
            c4 = c['ra_corn4']
            if c0  == 0:
                if c1 > 270:
                    c1 -= 360
                if c2 > 270:
                    c2 -= 360
                if c3 > 270:
                    c3 -= 360
                if c4 > 270:
                    c4 -= 360
            d1 = c['dec_corn1']
            d2 = c['dec_corn2']
            d3 = c['dec_corn3']
            d4 = c['dec_corn4']
            if distortion:
                c1 = c0 + (c1-c0) * cosfac
                c2 = c0 + (c2-c0) * cosfac
                c3 = c0 + (c3-c0) * cosfac
                c4 = c0 + (c4-c0) * cosfac
            ax.plot([c1, c2, c3, c4, c1],
                    [d1, d2, d3, d4, d1],color='skyblue')
        amin, amax = tile['ra_min'], tile['ra_max']
        dmin, dmax = tile['dec_min'], tile['dec_max']
        if c0 == 0:
            print('c0 is zero')
            if amin > 270:
                amin -= 360
            if amax > 270:
                amax -= 360
        r = np.arange(0,1.1,.1)
        r1 = np.ones(len(r))
        a1, d1 = r1 * amin, r * (dmax-dmin) + dmin
        a2, d2 = r * (amax-amin) + amin, r1 * dmax
        a3, d3 = r1 * amax, d1[::-1]
        a4, d4 = a2[::-1], r1*dmin
        aa = np.concatenate([a1,a2,a3,a4])
        dd = np.concatenate([d1,d2,d3,d4])
        #aa = [tile['ra_min'],tile['ra_min'],tile['ra_max'],tile['ra_max'],tile['ra_min']]
        #dd = [tile['dec_min'],tile['dec_max'],tile['dec_max'],tile['dec_min'],tile['dec_min']]
        if distortion:
            aa = c0 + (np.array(aa) - c0) * np.cos(dd * np.pi/180)
        ax.plot(aa,dd,color='red')
    ax.grid()
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    plt.show()
    return 1
