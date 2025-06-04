from hera import toolkitHome
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
from datetime import datetime, timedelta
import sklearn.metrics

tk = toolkitHome.getToolkit(toolkitName=toolkitHome.METEOROLOGY_LOWFREQ)

dirfile = r'/data3/Campaigns_Data/hera-ims/data/'
filelist = r'/data3/Campaigns_Data/hera-ims/data.json'
with open(filelist, 'r') as file:
    meta = json.load(file)

field='WD'
# field='WS'
# field='TD'
# field='RH'

tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY,projectName="topotesta")
# tk = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_RASTER_TOPOGRAPHY)

stations = []
for station in meta['deviceTypes'][0]['devices']:
    name = station['name']
    location = station['attributes'][2]['value']
    lat = location['latitude']
    lon = location['longitude']   
    elevation = tk.getPointElevation(lat,lon)
    try:
        data = pd.read_parquet(dirfile+name+'.parquet')
        # if ('WS' in data.keys()) and ('WD' in data.keys()):
        if (field in data.keys()):
           stations.append({'name':name,'lat':lat,'lon':lon,'elevation':elevation})
           print('-->', name)
        else:
           print('<', name)           
    except:
           print('-', name)           
      

best0=-9999
best1=0
best2=0        
scoretreshold = [0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
matcorr=np.zeros([len(stations),len(stations)])
for i in range(len(stations)):
    lat1=stations[i]['lat']*1000
    lon1=stations[i]['lon']*1000
    data1 = pd.read_parquet(dirfile+stations[i]['name']+'.parquet')
    data1=data1[data1[field]>-888]
    count=np.zeros(len(scoretreshold))
    for j in range(len(stations)):
        lat2=stations[j]['lat']*1000
        lon2=stations[j]['lon']*1000
        distscore=1
        distscore=((lat1-lat2)**2+(lon1-lon2)**2)/1000000.
        
        data2 = pd.read_parquet(dirfile+stations[j]['name']+'.parquet')
        data2=data2[data2[field]>-888]
        
        dfmerged = data1.merge(data2, on='datetime')
        
        maxspeed=2 #4 #2
        # dfmerged=dfmerged[dfmerged['WS_x']>maxspeed]
        if field=='TD':
            dfmerged[field+'_x']+=(stations[i]['elevation']-stations[j]['elevation'])/100
        
        if field=='WD':
            dfmerged.loc[(dfmerged[field+'_y'] - dfmerged[field+'_x'])>180,field+'_x']+=360
            dfmerged.loc[(dfmerged[field+'_x'] - dfmerged[field+'_y'])>180,field+'_y']+=360
            dfmerged.loc[(dfmerged[field+'_x'] + dfmerged[field+'_y'])<180,[field+'_x',field+'_y']]+=360
            # dfmerged=dfmerged[dfmerged['WS_x']>2]
        if len(dfmerged)>3600:
            # score = sklearn.metrics.r2_score(dfmerged[field+'_x'],dfmerged[field+'_y'])
            score = np.ma.corrcoef(np.ma.masked_invalid(dfmerged[field+'_x']), np.ma.masked_invalid(dfmerged[field+'_y']))[0].data[1]**1
            matcorr[i,j]=score*distscore
            if score>scoretreshold[5]:
                count+=1
                if score>best0 and i!=j:
                    best0=score
                    best1=i
                    best2=j
    stations[i]['corr']=np.sum(matcorr[i,:]>0.5)
    stations[i]['corr']=np.sum(matcorr[i,:]) # for distscore

for i in range(len(stations)):
    for j in range(len(scoretreshold)):
        stations[i]['corr'+str(j)]=np.sum(matcorr[i,:]>scoretreshold[j])
        # stations[i]['corr'+str(j)]=np.sum((matcorr.sum(axis=0)>scoretreshold[j]))
print(i,stations[i]['name'], count)
matcorr[matcorr==1]=np.nan

plt.figure()
plt.imshow(matcorr, origin='lower')
plt.colorbar()
plt.show()

# Extract values
lats = [d["lat"] for d in stations]
longs = [d["lon"] for d in stations]
corrs = [d["corr"] for d in stations]

# Create scatter plot
plt.figure(figsize=(8, 6))
sc = plt.scatter(longs, lats, c=corrs, cmap='coolwarm', edgecolors='black', s=100)
plt.colorbar(sc, label="corr5 Value")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Corr5 Values at Coordinates")

# Annotate each point with its corr5 value
for i, txt in enumerate(corrs):
    plt.text(longs[i], lats[i], f"{txt:.0f}", fontsize=9, ha='right')
# ax.set_aspect('equal', adjustable='datalim')
plt.axis('equal')  # This line sets the aspect ratio to 'equal'
plt.show()
        
##################################

def idw(x,y,stations, cuthi=False, s=0.5, alpha=2, elev=None, cam2=45, corr=False):
    for i in range(len(stations)):
        if stations[i][0] == x and stations[i][1] == y:
            return stations[i][2]

    wtotal=0
    for k in range(len(stations)):
        r = ((stations[k,0]-x)**2+(stations[k,1]-y)**2)  # it's r^2
        if r!=0:
            w=1./r**(alpha/2)
            wcuthi = 1.
            if cuthi:
                for m in range(len(stations)):
                    adist2 = float((stations[k][0] - stations[m][0]) ** 2 + (stations[k][1] - stations[m][1]) ** 2)
                    bdist2 = float((stations[m][0] - x) ** 2 + (stations[m][1] - y) ** 2)
                    if adist2>0 and bdist2>0:
                        cosgama = (adist2 + bdist2 - r) / (2.0 * adist2**.5 * bdist2**.5)
                        nangle = ((cosgama+1.)/2.)
                        if nangle<0:
                            if nangle<-0.01:
                                print ('nangle ngnalge')
                            else:
                                nangle = 0.
                                
                        wcuthi*=nangle**s
                w = w*wcuthi        
            if elev is not None:
                # cam2=45
                dz=elev-stations[k][3]
                welev=1./(1+cam2*dz**2)
                w*=welev
            if corr:
                # if stations[k,4]<=5:
                    # w=0
                w*=(1+stations[k,4]**0.5) # 0.5
                # w*=math.exp(-stations[k,4])
            wtotal+=w
    res=0
    for k in range(len(stations)):
        r = ((stations[k,0]-x)**2+(stations[k,1]-y)**2)
        if r == 0:
            res = stations[k,2]
        else:
            w=1./r**(alpha/2)
            wcuthi=1.
            if cuthi:
                for m in range(len(stations)):
                    adist2 = float((stations[k][0] - stations[m][0]) ** 2 + (stations[k][1] - stations[m][1]) ** 2)
                    bdist2 = float((stations[m][0] - x) ** 2 + (stations[m][1] - y) ** 2)
                    if adist2>0 and bdist2>0:
                        cosgama = (adist2 + bdist2 - r) / (2.0 * adist2**.5 * bdist2**.5)
                        nangle = ((cosgama+1.)/2.)
                        if nangle<0:
                            if nangle<-0.01:
                                print ('nangle ngnalge')
                            else:
                                nangle = 0.
                        wcuthi*=nangle**s
                w = w*wcuthi
            if elev is not None:
                # cam2=45
                dz=elev-stations[k][3]
                welev=1./(1+cam2*dz**2)
                w*=welev
                # print(elev,w,welev,cam2,dz)
            if corr:
                w*=(1+stations[k,4]**0.5) # 0.5
                # w*=math.exp(-stations[k,4])
                # if stations[k,4]<=5:
                    # w=0
            res+=stations[k,2]*w/wtotal
            
    return res

def crossvalidation(stations, cuthi = False, alpha = 2, s = None, elevation=False, cam2=45, corr=False, cyclic=False):
    
    r2 = 0
    r = 0
    rmsd = 0
    mae = 0
    size=len(stations)
    result = []    
    for i in range(size):
            # if i%100 == 0: print(i,'/',size)
            cv = []
            for j in range(size):
                if i!=j:
                   cv.append(stations[j])                      
            cv = np.asarray(cv)
            if elevation is False:
                result.append(idw(stations[i,0],stations[i,1],cv, cuthi=cuthi, alpha=alpha,corr=corr))
            else:
                result.append(idw(stations[i,0],stations[i,1],cv, cuthi=cuthi, alpha=alpha, elev=stations[i,3], cam2=cam2,corr=corr))
                   
    result=np.asarray(result)
    if cyclic:
        result[result>360]=result[result>360]-360
        result[result-stations[:,2]>180] = result[result-stations[:,2]>180] - 360
        result[stations[:,2]-result>180] = result[stations[:,2]-result>180] + 360
                
    r2   += sklearn.metrics.r2_score(stations[:,2], result)
    rmsd += np.sqrt(((stations[:,2] - result)   ** 2).mean()) 
    mae += (np.abs(stations[:,2] - result).mean()) 
    r   += np.ma.corrcoef(np.ma.masked_invalid(stations[:,2]), np.ma.masked_invalid(result))[0].data[1]
        
    return round(r2,3), round(rmsd,3), round(mae,3) , round(r,3)

def totalcrossvalidation(stations, cuthi = False, alpha = 2, corr=False):
    
    result = []    
    for i in range(len(stations)):
            # if i%100 == 0: print(i,'/',size)
            cv = []
            for j in range(len(stations)):
                if i!=j:
                   cv.append(stations[j])                      
            cv = np.asarray(cv)
            result.append(idw(stations[i,0],stations[i,1],cv, cuthi=cuthi, alpha=alpha,corr=corr))
                                           
    return result

kind='WD'
col='corr'
for i in range(len(stations)):
    if stations[i]['lat']<100:
        stations[i]['lat']*=1000
        stations[i]['lon']*=1000
    print(i, stations[i]['lat'], stations[i]['lon'])

ranks0=np.zeros(2)
ranks1=np.zeros(2)
ranks2=np.zeros(2)
ranks3=np.zeros(2)
ranksrmsd=[]
ranksr2=[]
ranksr=[]
ranksmae=[]
time_string = '2021-04-01T00:00:00+03:00'
time_obj = datetime.fromisoformat(time_string)
for t in range(12):
    stwss=[]
    # stwds=[]
    for i in range(len(stations)):
        # print(i, len(stations))
        try:
            data1 = pd.read_parquet(dirfile+stations[i]['name']+'.parquet')
            data1=data1[data1[field]>-888]
            new_time_obj = time_obj + timedelta(hours=t*3)
            new_time_string = new_time_obj.isoformat()        
            stws=data1[data1['datetime']==new_time_string][[field]].values[0][0]
            if field=='TD':
                stws+=stations[i]['elevation']/100
            # if data1[data1['datetime']==new_time_string][['WS']].values[0][0]>2:
            stwss.append([stations[i]['lat'],stations[i]['lon'],stws,stations[i]['elevation'],stations[i][col]])
            # stwds.append([stations[i]['lat'],stations[i]['lon'],stwd,stations[i]['elevation'],stations[i]['corr']])
        except:
            continue
    # print(len(stwss))
    stwss = np.asarray(stwss)
    # stwds = np.asarray(stwds)
    # stwds = numpy.vstack([stwds, [33,35,270, 333,1]])

    # print(t,'**************'+new_time_string, ranks)
    if kind=='WD':
        cyc=True
    else:
        cyc=False
        
    rank=[]
    rank.append(crossvalidation(stwss, cyclic=cyc))
    rank.append(crossvalidation(stwss,corr=True, cyclic=cyc))
    # rank.append(crossvalidation(stwss,cuthi=True, cyclic=cyc))
    # rank.append(crossvalidation(stwss,cuthi=True,corr=True, cyclic=cyc))
    # rank.append(crossvalidation(stwss, elevation=True,cam2=50/10000, cuthi=True)) # 5/10000
    # rank.append(crossvalidation(stwss, elevation=True,cam2=50/10000, cuthi=True, corr=True)) # 5/10000
    # rank.append(crossvalidation(stwss, elevation=True,cam2=50/10000)) # 5/10000
    # rank.append(crossvalidation(stwss, elevation=True,cam2=50/10000,corr=True)) # 5/10000
    rank=np.asarray(rank)
    ranks0[rank[:,0].argmax()]+=1 # r2
    ranks1[rank[:,1].argmin()]+=1 # rmsd
    ranks2[rank[:,2].argmin()]+=1 # mae
    ranks3[rank[:,3].argmax()]+=1 # r
    ranksmae.append(rank[:,2])
    ranksrmsd.append(rank[:,1])
    ranksr2.append(rank[:,0])
    ranksr.append(rank[:,3])
    # print(rank[:,1].argmin(),rank)
    rmsdtemp = np.asarray(ranksrmsd)
# rmsd, r2, mae, r
    print(t,'*******'+new_time_string, ranks0,rmsdtemp.sum(axis=0),np.asarray(ranksr2).sum(axis=0),np.asarray(ranksmae).sum(axis=0),np.asarray(ranksr).sum(axis=0))
    

time_string = '2021-04-01T00:00:00+03:00'
time_obj = datetime.fromisoformat(time_string)
rsidw=[]
rscorr=[]
rscuthi=[]
rscuthicorr=[]
obs=[]
for t in range(188887):
    stwss=[]
    # stwds=[]
    for i in range(len(stations)):
        # print(i, len(stations))
        try:
            data1 = pd.read_parquet(dirfile+stations[i]['name']+'.parquet')
            data1=data1[data1[field]>-888]
            new_time_obj = time_obj + timedelta(hours=t*1)
            new_time_string = new_time_obj.isoformat()        
            stws=data1[data1['datetime']==new_time_string][[field]].values[0][0]
            if field=='TD':
                stws+=stations[i]['elevation']/100
            # if data1[data1['datetime']==new_time_string][['WS']].values[0][0]>2:
            stwss.append([stations[i]['lat'],stations[i]['lon'],stws,stations[i]['elevation'],stations[i][col]])
            # stwds.append([stations[i]['lat'],stations[i]['lon'],stwd,stations[i]['elevation'],stations[i]['corr']])
        except:
            continue
    # print(len(stwss))
    stwss = np.asarray(stwss)
    stwss[:,4][stwss[:,4]<-100.0]=-100.0
    stwss[:,4]-=stwss[:,4].min()
    rsidw += totalcrossvalidation(stwss)
    rscorr += totalcrossvalidation(stwss,corr=True)
    rscuthi += totalcrossvalidation(stwss, cuthi=True)
    rscuthicorr += totalcrossvalidation(stwss,cuthi=True,corr=True)
    obs += stwss[:,2].tolist()
    if t % 10 == 0:
        print(t, new_time_obj)
    if t % 10 == 0:
        rsidw1=np.asarray(rsidw)    
        rscorr1=np.asarray(rscorr)    
        rscuthi1=np.asarray(rscuthi)    
        rscuthicorr1=np.asarray(rscuthicorr)    
        obs1=np.asarray(obs)    
        # print('r2, idw, corr:', round(sklearn.metrics.r2_score(obs1, rsidw1),2), round(sklearn.metrics.r2_score(obs1, rscorr1),2))
        # print('rmsd, idw, corr:', round(np.sqrt(((obs1 - rsidw1) ** 2).mean()),2), round(np.sqrt(((obs1 - rscorr1)   ** 2).mean()),2))
        # print('mae, idw, corr:', round(np.abs(obs1 - rsidw1).mean(),2), round(np.abs(obs1 - rscorr1).mean(),2))
        # print('bias, idw, corr:', round((obs1 - rsidw1).mean(),2), round((obs1 - rscorr1).mean(),2))
        # print('r, idw, corr:', round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rsidw1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscorr1))[0].data[1],2))
        print('r2, idw, corr, cuthi, cuthicorr:', round(sklearn.metrics.r2_score(obs1, rsidw1),2), round(sklearn.metrics.r2_score(obs1, rscorr1),2), round(sklearn.metrics.r2_score(obs1, rscuthi1),2), round(sklearn.metrics.r2_score(obs1, rscuthicorr1),2))
        print('rmsd, idw, corr, cuthi, cuthicorr:', round(np.sqrt(((obs1 - rsidw1) ** 2).mean()),2), round(np.sqrt(((obs1 - rscorr1)   ** 2).mean()),2), round(np.sqrt(((obs1 - rscuthi1) ** 2).mean()),2), round(np.sqrt(((obs1 - rscuthicorr1)   ** 2).mean()),2))
        print('mae, idw, corr, cuthi, cuthicorr:', round(np.abs(obs1 - rsidw1).mean(),2), round(np.abs(obs1 - rscorr1).mean(),2), round(np.abs(obs1 - rscuthi1).mean(),2), round(np.abs(obs1 - rscuthicorr1).mean(),2))
        print('bias, idw, corr, cuthi, cuthicorr:', round((obs1 - rsidw1).mean(),2), round((obs1 - rscorr1).mean(),2), round((obs1 - rscuthi1).mean(),2), round((obs1 - rscuthicorr1).mean(),2))
        print('r, idw, corr, cuthi, cuthicorr:', round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rsidw1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscorr1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscuthi1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscuthicorr1))[0].data[1],2))
        if kind=='WD':
            print('cyc')
            rsidw1[rsidw1>360]=rsidw1[rsidw1>360]-360
            rsidw1[rsidw1-obs1>180] = rsidw1[rsidw1-obs1>180] - 360
            rsidw1[obs1-rsidw1>180] = rsidw1[obs1-rsidw1>180] + 360    
            rscorr1[rscorr1>360]=rscorr1[rscorr1>360]-360
            rscorr1[rscorr1-obs1>180] = rscorr1[rscorr1-obs1>180] - 360
            rscorr1[obs1-rscorr1>180] = rscorr1[obs1-rscorr1>180] + 360

            rscuthi1[rscuthi1>360]=rscuthi1[rscuthi1>360]-360
            rscuthi1[rscuthi1-obs1>180] = rscuthi1[rscuthi1-obs1>180] - 360
            rscuthi1[obs1-rscuthi1>180] = rscuthi1[obs1-rscuthi1>180] + 360
            rscuthicorr1[rscuthicorr1>360]=rscuthicorr1[rscuthicorr1>360]-360
            rscuthicorr1[rscuthicorr1-obs1>180] = rscuthicorr1[rscuthicorr1-obs1>180] - 360
            rscuthicorr1[obs1-rscuthicorr1>180] = rscuthicorr1[obs1-rscuthicorr1>180] + 360

            print('r2, idw, corr, cuthi, cuthicorr:', round(sklearn.metrics.r2_score(obs1, rsidw1),2), round(sklearn.metrics.r2_score(obs1, rscorr1),2), round(sklearn.metrics.r2_score(obs1, rscuthi1),2), round(sklearn.metrics.r2_score(obs1, rscuthicorr1),2))
            print('rmsd, idw, corr, cuthi, cuthicorr:', round(np.sqrt(((obs1 - rsidw1) ** 2).mean()),2), round(np.sqrt(((obs1 - rscorr1)   ** 2).mean()),2), round(np.sqrt(((obs1 - rscuthi1) ** 2).mean()),2), round(np.sqrt(((obs1 - rscuthicorr1)   ** 2).mean()),2))
            print('mae, idw, corr, cuthi, cuthicorr:', round(np.abs(obs1 - rsidw1).mean(),2), round(np.abs(obs1 - rscorr1).mean(),2), round(np.abs(obs1 - rscuthi1).mean(),2), round(np.abs(obs1 - rscuthicorr1).mean(),2))
            print('bias, idw, corr, cuthi, cuthicorr:', round((obs1 - rsidw1).mean(),2), round((obs1 - rscorr1).mean(),2), round((obs1 - rscuthi1).mean(),2), round((obs1 - rscuthicorr1).mean(),2))
            print('r, idw, corr, cuthi, cuthicorr:', round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rsidw1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscorr1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscuthi1))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs1), np.ma.masked_invalid(rscuthicorr1))[0].data[1],2))
        
    
rsidw=np.asarray(rsidw)    
rscorr=np.asarray(rscorr)    
obs=np.asarray(obs)    
print('r2, idw, corr:', round(sklearn.metrics.r2_score(obs, rsidw),2), round(sklearn.metrics.r2_score(obs, rscorr),2))
print('rmsd, idw, corr:', round(np.sqrt(((obs - rsidw) ** 2).mean()),2), round(np.sqrt(((obs - rscorr)   ** 2).mean()),2))
print('mae, idw, corr:', round(np.abs(obs - rsidw).mean(),2), round(np.abs(obs - rscorr).mean(),2))
print('bias, idw, corr:', round((obs - rsidw).mean(),2), round((obs - rscorr).mean(),2))
print('r, idw, corr:', round(np.ma.corrcoef(np.ma.masked_invalid(obs), np.ma.masked_invalid(rsidw))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs), np.ma.masked_invalid(rscorr))[0].data[1],2))
if kind=='WD':
    print('cyc')
    rsidw[rsidw>360]=rsidw[rsidw>360]-360
    rsidw[rsidw-obs>180] = rsidw[rsidw-obs>180] - 360
    rsidw[obs-rsidw>180] = rsidw[obs-rsidw>180] + 360    
    rscorr[rscorr>360]=rscorr[rscorr>360]-360
    rscorr[rscorr-obs>180] = rscorr[rscorr-obs>180] - 360
    rscorr[obs-rscorr>180] = rscorr[obs-rscorr>180] + 360
    print('r2, idw, corr:', round(sklearn.metrics.r2_score(obs, rsidw),2), round(sklearn.metrics.r2_score(obs, rscorr),2))
    print('rmsd, idw, corr:', round(np.sqrt(((obs - rsidw) ** 2).mean()),2), round(np.sqrt(((obs - rscorr)   ** 2).mean()),2))
    print('mae, idw, corr:', round(np.abs(obs - rsidw).mean(),2), round(np.abs(obs - rscorr).mean(),2))
    print('bias, idw, corr:', round((obs - rsidw).mean(),2), round((obs - rscorr).mean(),2))
    print('r, idw, corr:', round(np.ma.corrcoef(np.ma.masked_invalid(obs), np.ma.masked_invalid(rsidw))[0].data[1],2), round(np.ma.corrcoef(np.ma.masked_invalid(obs), np.ma.masked_invalid(rscorr))[0].data[1],2))


plt.figure()
plt.scatter(obs[::100],rsidw[::100],s=0.1,c='r')
plt.title('idw')
plt.show()

plt.figure()
plt.scatter(obs[::10],rscorr[::10],s=0.5,c='b')
plt.scatter(obs[::10],rsidw[::10],s=0.5,c='r')
plt.title('corr')
plt.show()

plt.figure()
plt.scatter(rsidw[::10],rscorr[::10],s=0.1)
plt.title('corr')
plt.show()

ranksrmsd = np.asarray(ranksrmsd)
print(ranksrmsd.sum(axis=0))
print('if correlation is zero change it to one!!!!')

plt.figure()
plt.plot(stwss[:,2],stwss[:,4],'*')
plt.title(col)
plt.show()

# crossvalidation(stwss, elevation=True,cam2=9/10000) # 5/10000
# crossvalidation(stwss, elevation=True,cam2=60/10000) # 5/10000
