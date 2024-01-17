def







def kmeans_v3(v,K):
    data_set = np.dstack(v.flatten())
    model = KMeans(n_clusters=K)
    model.fit(data_set[0].T)
    return model.labels_

seg_tools.init()
I = mpimg.imread('implant.bmp')
dim_y, dim_x = np.shape(I)
v = np.reshape(I.astype(float),(dim_y*dim_x,1))
c = np.array([50,150,200])
t1 = time.perf_counter()
for n in range(0,10):
    label1, i =kmeans_v1(v,c)
print('V1 : iterations :',i)
t2 = time.perf_counter()
print('V1 : temps :',(t2-t1)/10)
for n in range(0,1000):
    label2, i =kmeans_v1(v,c)
print('V2:iterations:',i)
