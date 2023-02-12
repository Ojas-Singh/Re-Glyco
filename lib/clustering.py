import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lib import graph
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stat
from pylab import cm

def normalizetorsion(df):
    tor = df.loc[:, df.columns!='i'].to_numpy()
    tor=tor.T
    normalized = []
    
    for angle in range(len(tor)): 
        x = 0
        y = 0
        for i in tor[angle]:
            x += np.cos(np.deg2rad(i))
            y += np.sin(np.deg2rad(i))
        average_angle = np.arctan2(y, x)
        # print(average_angle)
        for i in range(len(tor[angle])):
            tor[angle][i]=np.arctan(np.tan(np.deg2rad(tor[angle][i])-average_angle))
    tor=tor.T
    return tor

def pcawithT(tor):
    
    pca = PCA(n_components=2)
    t = pca.fit_transform(tor)
    x=[]
    y=[]
    for i in t:
        x.append(i[0])
        y.append(i[1])
    df = pd.DataFrame(data=np.column_stack((x,y)), columns = ['X','Y'])
    PCA_components = pd.DataFrame(t)
    return df,PCA_components


def pcawithG(frames,idx_noH):
    G = np.zeros((len(frames),int(len(frames[0][np.asarray(idx_noH,dtype=int)])*(len(frames[0][np.asarray(idx_noH,dtype=int)])+1)/2)))
    for i in range(len(frames)):
        G[i]= graph.G_flatten(frames[i][np.asarray(idx_noH,dtype=int)])
    pca = PCA(n_components=20)
    t = pca.fit_transform(G)
    x=[]
    y=[]
    for i in t:
        x.append(i[0])
        y.append(i[1])
    df = pd.DataFrame(data=np.column_stack((x,y)), columns = ['X','Y'])
    
    PCA_components = pd.DataFrame(t)
    return df,PCA_components


def findmaxima(f):
    f = -1*f
    a = np.pad(f, (1, 1), mode='constant',
            constant_values=(np.amax(f), np.amax(f)))
    loc_min = []
    rows = a.shape[0]
    cols = a.shape[1]
    for ix in range(0, rows - 1):
        for iy in range(0, cols - 1):
                    if (a[ix, iy] < a[ix, iy + 1]
                        and a[ix, iy] < a[ix, iy - 1]
                        and a[ix, iy] < a[ix + 1, iy]
                        and a[ix, iy] < a[ix + 1, iy - 1]
                        and a[ix, iy] < a[ix + 1, iy + 1]
                        and a[ix, iy] < a[ix - 1, iy]
                        and a[ix, iy] < a[ix - 1, iy - 1]
                        and a[ix, iy] < a[ix - 1, iy + 1]):
                        temp_pos = (ix-1, iy-1)
                        loc_min.append(temp_pos)
    return loc_min

def nux(x,d,xmin):
    # return xmin + d*x/100
    return (x-xmin)*100/d

def nuy(y,d,ymin):
    # return ymin + d*y/100
    return (y-ymin)*100/d

def filterlow(data,nth):
    x=[]
    y=[]
    for i in range(len(data)):
        if i%nth==0:
            x.append(data.iloc[i,0])
            y.append(data.iloc[i,1])
    s=pd.DataFrame([x,y])
    s=s.transpose()
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = stat.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    fig = plt.figure()
    mpl.rcParams['font.family'] = 'Cambria'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 2
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(xx, yy, f, cmap=plt.cm.YlGnBu_r)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title("PCA Space Density")
    ax.view_init(elev=-15, azim=-59)
    ax.contourf(xx, yy, f, zdir='z', offset=-0.8, cmap=plt.cm.YlGnBu_r)
    # plt.show()
    plt.savefig('PCA_KDE.png',dpi=450)
    plt.clf()
    l=[]
    for i in range(len(x)):
        l.append([f[int(nux(x[i],np.abs(xmax-xmin),xmin))-1][int(nuy(y[i],np.abs(ymax-ymin),ymin))-1],i])
    l.sort()
    
    
    idx_top=np.ones(len(x),dtype=bool) 
    idx_bottom = np.zeros(len(x),dtype=bool) 
    for i in range(int(.1*len(x))):
        idx_top[l[i][1]] = False
        idx_bottom[l[i][1]] = True
    x= np.asarray(x)
    y= np.asarray(y)
    fig = plt.figure()
    ax = fig.gca()
    cfset = ax.contourf(xx, yy, f, cmap='Blues')
    ax.scatter(x[idx_top],y[idx_top],color="#78517C",s=.2)
    ax.scatter(x[idx_bottom],y[idx_bottom],color="#F65058FF",s=.2)
    ax.set_title("Conformation Filter (>10%)")
    plt.savefig('PCA_filter.png',dpi=450)
    
    return idx_top

def cluster(data,numofcluster,idx):
    from sklearn.cluster import KMeans,SpectralCoclustering
    x=[]
    y=[]
    for i in range(len(data)):
        if i%1==0:
            x.append(data.iloc[i,0])
            y.append(data.iloc[i,1])
    s=pd.DataFrame([x[idx],y[idx]])
    s=s.transpose()
    ini=[]
    for i in range(numofcluster):
        ini.append([xw[i],yw[i]])
    clustering = KMeans(n_clusters=numofcluster,init=ini,n_init=1,max_iter=1).fit(s)
    clustering_labels= clustering.labels_
    ax.scatter(s.iloc[:,0],s.iloc[:,1],s=0.1,alpha=0.4,c=clustering_labels.astype(float))
    ax.scatter(xw,yw,color="blue")
    plt.show()

