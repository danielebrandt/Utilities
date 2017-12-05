import numpy as np
from numpy import random

def inc2latVec(Vinc,G2,G3):
	Vlat = np.zeros(len(Vinc))
	for i in range(len(Vinc)):
		Vlat[i] = inc2lat(Vinc[i],G2,G3)
	return Vlat 
def lat2incVec(Vlat, G2, G3):
	Vinc = np.zeros(len(Vlat))
	for i in range(len(Vlat)):
		Vinc[i] = lat2inc(Vlat[i],G2,G3)
	return Vinc
def inc2lat(inc, G2, G3):
    """Determinates the latitudes values from a magnetic inclination 
    inc is only one value of inclination (float)
	G2 can varies: 0, 0.1, 0.2 ...
    G3 can be less then 0.3, for higher values there are more than one solution, so it is not possible   
	From tan(INC) = A(teta)/B(teta) Merril et al. (1996) pag. 232 we find teta (colatitude) finding the root from H = A-Btan(Inc)
    using the  bisection method, also known as  dichotomy method"""
    
    parar = 0
    teta1 = 0.
    teta2 = 180.

    while (parar!=1):
        tetam = (teta1+teta2)/2
        
        if (abs(fun_h(tetam,inc, G2,G3))<0.00000001 or (tetam-teta1)<0.00000001) :
            parar = 1
        if (fun_h(tetam,inc, G2,G3)*fun_h(teta1,inc, G2,G3)>0):
            teta1=tetam
        if (fun_h(tetam,inc, G2,G3)*fun_h(teta2,inc, G2,G3)>0):
            teta2=tetam

    lat = 90 - tetam
    return lat
def lat2inc(lat,G2,G3):
    """Diretc application of the function tan(INC) = A(teta)/B(teta) for finding an value of inclination from a latitude value
    we can use any value of contribution of G2 e G3: 0, 0.1, 0.2, 0.3, 0.4...
    Merril et al. (1996) pag. 232"""
    teta = 90. - lat
    r_te = np.deg2rad(teta) #radianos do teta
    costeta = np.cos(r_te)
    cos2teta = costeta*costeta
    cos3teta = cos2teta*costeta
    
    sinteta = np.sin(r_te)

    A = 2*costeta + G2*(4.5*cos2teta - 1.5) + G3*(10*cos3teta - 6*costeta)
    B = sinteta + G2*(3*costeta*sinteta) + G3*(7.5*cos2teta*sinteta - 1.5*sinteta)
    inc = np.rad2deg(np.arctan(A/B))


    return inc
def fun_h(teta, inc, G2, G3):
    """ Finding the root of the function h(teta) we find the teta that corresponds a one value os inclination "inc" 
    H = A-Btan(Inc) from tan(INC) = A(teta)/B(teta) Merril et al. (1996) pag. 232"""
    r_te = np.deg2rad(teta) #radianos do teta
    costeta = np.cos(r_te)
    cos2teta = costeta*costeta
    cos3teta = cos2teta*costeta
    cos4teta = cos3teta*costeta
    sinteta = np.sin(r_te)
    sin2teta = sinteta*sinteta
    tanI = np.tan(np.deg2rad(inc)) 

    A = 2*costeta + G2*(4.5*cos2teta - 1.5) + G3*(10*cos3teta - 6*costeta)
    B = sinteta + G2*(3*costeta*sinteta) + G3*(7.5*cos2teta*sinteta - 1.5*sinteta) 

    h = A - B*tanI
    
    return h

def dir2vgp_oneplace(DI, Lats, Lons):
	"""Receive one matrix with 2 columms (Dec, Inc) and the values Lats, Lons of one place. Gives a matrix with two columms (PGVLon, PGVlat)
	without alfa95, dp and dm"""
	LonLat = np.zeros((len(DI),2)) #crio uma matriz com 2 colunas e comprimento N 
	for i in range(0,len(DI)):
		Dec_r = np.deg2rad(DI[i,0])
		Inc_r = np.deg2rad(DI[i,1])
        
		Lats_r = np.deg2rad(Lats)
		Lons_r = np.deg2rad(Lons)        
        
		paleolat_r = np.arctan2(np.tan(Inc_r), 2.)
		p_r = (np.pi/2.) - paleolat_r
		Latp = np.rad2deg(np.arcsin(np.sin(Lats_r) * np.cos(p_r) + (np.cos(Lats_r) * np.sin(p_r) * np.cos(Dec_r))))
		Beta =  (np.sin(p_r) * np.sin(Dec_r)) / np.cos(np.deg2rad(Latp))
		if Beta>1.:
			Beta = 1.
		if Beta<-1.:
			Beta = -1.
		Beta = np.arcsin(Beta)
		if  np.isnan(Beta)== True:
			print ('Problema com o Beta = % ' % Beta)
        
		if np.cos(p_r) >= (np.sin(np.deg2rad(Latp)) * np.sin(Lats_r)):
			Lonp = np.rad2deg(Lons_r + Beta)
		if np.cos(p_r) < (np.sin(np.deg2rad(Latp)) * np.sin(Lats_r)):
			Lonp = np.rad2deg(Lons_r + np.pi - Beta)
		if Lonp > 360: Lonp = Lonp - 360.
		if Lonp < 0: Lonp = 360. + Lonp
		if np.isnan(Lonp)== True:
			print ('Problema com a Lonp = % ' % Lonp)
		LonLat[i,0] = Lonp
		if np.isnan(Latp)== True:
			print ('Problema com a Lonp = % ' % Latp)
		LonLat[i,1] = Latp 
        
	return LonLat
def dir2vgp_severalplaces(DI_5):
	"""Generates a matriz de with N lines and 4 columns Plon, Plat, dm, dp from a set of data Dec,Inc,alfa,Lats,Lons"""
	LonLatdmdp = np.zeros((len(DI_5),4)) #crio uma matriz com 4 colunas e comprimento N 
	for i in range(0,len(DI_5)):
		Dec_r = np.deg2rad(DI_5[i,0])
		Inc_r = np.deg2rad(DI_5[i,1])
		alfa = DI_5[i,2]
		Lats_r = np.deg2rad(DI_5[i,3])
		Lons_r = np.deg2rad(DI_5[i,4])        
		paleolat_r = np.arctan2(np.tan(Inc_r), 2.)
		p_r = (np.pi/2.) - paleolat_r
		Latp = np.rad2deg(np.arcsin(np.sin(Lats_r) * np.cos(p_r) + (np.cos(Lats_r) * np.sin(p_r) * np.cos(Dec_r))))
		Beta =  (np.sin(p_r) * np.sin(Dec_r)) / np.cos(np.deg2rad(Latp))
		if Beta>1.:
			Beta = 1.
		if Beta<-1.:
			Beta = -1.
		Beta = np.arcsin(Beta)
		if  np.isnan(Beta)== True:
			print ('Problema com o Beta = % ' % Beta)
		if np.cos(p_r) >= (np.sin(np.deg2rad(Latp)) * np.sin(Lats_r)):
			Lonp = np.rad2deg(Lons_r + Beta)
		if np.cos(p_r) < (np.sin(np.deg2rad(Latp)) * np.sin(Lats_r)):
			Lonp = np.rad2deg(Lons_r + np.pi - Beta)
		if Lonp > 360: Lonp = Lonp - 360
		if Lonp < 0: Lonp = 360 + Lonp
		dm = alfa * np.sin(p_r) / np.cos(Inc_r)
		dp = 0.5 * alfa * (1 + (3 * (np.cos(p_r)**2)))
		LonLatdmdp[i,0] = Lonp
		LonLatdmdp[i,1] = Latp 
		LonLatdmdp[i,2] = dm
		LonLatdmdp[i,3] = dp  
	return LonLatdmdp
def dispS(LonLat):
	"""Calculates the dispersion S from a set of PGVs LonLat"""
	N=len(LonLat)#length lista de dados eh N
	LoLaXYZ = np.zeros((N,4))# XYZ and deltas i
	M= np.zeros((N,3))#matriz usada pra XYZ sem deltai
	for j in range(0,N):
		#coordenadas esf. para ortogonais
		LoLaXYZ[j,0] =  np.cos(np.deg2rad(LonLat[j,0]))*np.cos(np.deg2rad(LonLat[j,1]))
		LoLaXYZ[j,1] =  np.sin(np.deg2rad(LonLat[j,0]))*np.cos(np.deg2rad(LonLat[j,1]))
		LoLaXYZ[j,2] =  np.sin(np.deg2rad(LonLat[j,1]))
		LoLaXYZ[j,3] = 0.0
		#print "LoLaXYZ"
		#imprima_matriz(LoLaXYZ)
	M = LoLaXYZ[:,0:3]
	mediaXYZ = [sum(M[:,0]),sum(M[:,1]),sum(M[:,2])]
	#print "media", mediaXYZ[0],mediaXYZ[1], mediaXYZ[2]
	
	R = np.sqrt(mediaXYZ[0]*mediaXYZ[0]+mediaXYZ[1]*mediaXYZ[1]+mediaXYZ[2]*mediaXYZ[2])
	#print "R",R
	mediaXYZn =  mediaXYZ/R #Media de Fisher em XYZ
	LoLaXYZ[:,3] = np.rad2deg(np.arccos(np.matmul(M,mediaXYZn))) #deltais
	for j in range(0,N):
		if abs(LoLaXYZ[j,3])>90.:
			LoLaXYZ[j,3]=180.-abs(LoLaXYZ[j,3]) #for reversed data 
	S=np.sqrt(sum(LoLaXYZ[:,3]*LoLaXYZ[:,3])/(N-1))#calculo da dispersao S  
	return S

def dispS_boot(VGPS, NB):
    """Calculates de bootstrap 95 percent of confidence regions of dispersion S"""
    low = int(.025 * NB)
    high = int(.975 * NB)
    S_b = np.zeros(NB)
    for i in range(NB):
        VGPS_B = np.zeros([len(VGPS),2])
        for k in range(len(VGPS)):
            random.seed()
            ind = random.randint(len(VGPS))
            VGPS_B[k,:]=VGPS[ind,:]
        S_b[i] = dispS(VGPS_B)
    S_b.sort()
      
    return S_b[low],S_b[high]

def Scut(LonLat,boot=0,vand=1, A=0., NB = 1000):
	"""Calculates the dispersion S using a cutoff angle A (from Vandamme 1994), returns final number of VGPs, S and A"""
	N=len(LonLat)#size of data set
	if vand == 1:
		A = 0. #cutoff angle
	Deltamax = 180.#maximum Deltai 
	LoLaXYZ = np.zeros((N,6))# XYZ and deltas i Lon Lat em degrees
	Lola = np.zeros([N,2])
	M= np.zeros((N,3))#only XYZ without deltai
	for j in range(0,N):
		LoLaXYZ[j,0] =  np.cos(np.deg2rad(LonLat[j,0]))*np.cos(np.deg2rad(LonLat[j,1]))
		LoLaXYZ[j,1] =  np.sin(np.deg2rad(LonLat[j,0]))*np.cos(np.deg2rad(LonLat[j,1]))
		LoLaXYZ[j,2] =  np.sin(np.deg2rad(LonLat[j,1]))
		LoLaXYZ[j,3] = 0.0
		LoLaXYZ[j,4] = LonLat[j,0]
		LoLaXYZ[j,5] = LonLat[j,1]
	while Deltamax >A:
		M = LoLaXYZ[:,0:3]
		mediaXYZ = [sum(M[:,0]),sum(M[:,1]),sum(M[:,2])]
		R = np.sqrt(mediaXYZ[0]*mediaXYZ[0]+mediaXYZ[1]*mediaXYZ[1]+mediaXYZ[2]*mediaXYZ[2])
		mediaXYZn =  mediaXYZ/R #Media de Fisher em XYZ
		LoLaXYZ[:,3] = np.rad2deg(np.arccos(np.matmul(M,mediaXYZn))) #deltais
		for j in range(0,N):
			if LoLaXYZ[j,3]>90.:
				LoLaXYZ[j,3]=180.-LoLaXYZ[j,3] #for reversed data 
		seq = np.argsort(LoLaXYZ[:,3])#sorting using the deltai values
		LoLaXYZ = LoLaXYZ[seq]#sorting using the deltai values
		S=np.sqrt(sum(LoLaXYZ[:,3]*LoLaXYZ[:,3])/(N-1))# dispersion S  
		if vand==1:
			A=1.8*S+5 #cutoff angle 
		Deltamax = LoLaXYZ[N-1,3]#deltai maximum at the end of the list
		if Deltamax > A:
			N=N-1
			LoLaXYZ = LoLaXYZ[0:N] #removes the last line
	if boot==1:
		S_b_L,S_b_H = dispS_boot(LoLaXYZ[:,4:6], NB)
		return N, S, S_b_L,S_b_H, A
	else:
		return N, S, A


def dispS_boot(VGPS, NB):
    """Calculates de bootstrap 95 percent of confidence regions of dispersion S"""
    low = int(.025 * NB)
    high = int(.975 * NB)
    S_b = np.zeros(NB)
    for i in range(NB):
        VGPS_B = np.zeros([len(VGPS),2])
        for k in range(len(VGPS)):
            random.seed()
            ind = random.randint(len(VGPS))
            VGPS_B[k,:]=VGPS[ind,:]
        S_b[i] = dispS(VGPS_B)
    S_b.sort()
      
    return S_b[low],S_b[high]
def Cut_VGPS(LonLat, vand, A):
    """Removes VGPS using cutoff of Vandamme vand = 1, or fixed cutoff, returns an array with 
    the collums lon, lat and 0 or 1: 0 for cut data, 1 accepted data"""
    N=len(LonLat)#size of data set
    if vand == 1:
        A = 0.

    Deltamax = 180.#maximum Deltai 
    LoLaXYZ = np.zeros((N,7))# XYZ and deltas i Lon Lat em degrees
    Lola = np.zeros([N,2])
    M= np.zeros((N,3))#only XYZ without deltai
    for j in range(0,N):
        LoLaXYZ[j,0] =  np.cos(np.deg2rad(LonLat[j,0]))*np.cos(np.deg2rad(LonLat[j,1]))
        LoLaXYZ[j,1] =  np.sin(np.deg2rad(LonLat[j,0]))*np.cos(np.deg2rad(LonLat[j,1]))
        LoLaXYZ[j,2] =  np.sin(np.deg2rad(LonLat[j,1]))
        LoLaXYZ[j,3] = 0.0
        LoLaXYZ[j,4] = LonLat[j,0]
        LoLaXYZ[j,5] = LonLat[j,1]
        LoLaXYZ[j,6] = j
    while Deltamax >A:
        M = LoLaXYZ[:,0:3]
        mediaXYZ = [sum(M[:,0]),sum(M[:,1]),sum(M[:,2])]
        R = np.sqrt(mediaXYZ[0]*mediaXYZ[0]+mediaXYZ[1]*mediaXYZ[1]+mediaXYZ[2]*mediaXYZ[2])
        mediaXYZn =  mediaXYZ/R #Media de Fisher em XYZ
        LoLaXYZ[:,3] = np.rad2deg(np.arccos(np.matmul(M,mediaXYZn))) #deltais
        for j in range(0,N):
            if LoLaXYZ[j,3]>90.:
                LoLaXYZ[j,3]=180.-LoLaXYZ[j,3] #for reversed data 
        seq = np.argsort(LoLaXYZ[:,3])#sorting using the deltai values
        LoLaXYZ = LoLaXYZ[seq]#sorting using the deltai values
        S=np.sqrt(sum(LoLaXYZ[:,3]*LoLaXYZ[:,3])/(N-1))# dispersion S  
        if vand==1:
            A=1.8*S+5 #cutoff angle 
        Deltamax = LoLaXYZ[N-1,3]#deltai maximum at the end of the list
        if Deltamax > A:
            N=N-1
            LoLaXYZ = LoLaXYZ[0:N] #removes the last line
    LonLatv = np.zeros([len(LonLat),3])
    for j in range(len(LonLat)):
        for w in range(len(LoLaXYZ)):
            if LoLaXYZ[w,6]==j:
                LonLatv[j,0]=LoLaXYZ[w,4]
                LonLatv[j,1]=LoLaXYZ[w,5]
                LonLatv[j,2]=1
    return LonLatv


