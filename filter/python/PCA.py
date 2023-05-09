import numpy as np

class PCA():
    """
        input has dimensions [n_components,n_samples]
    """
    def __init__(self, input):
        print('Building PCA object...')
        inp = np.zeros(input.shape)
        mean = np.mean(input,axis=1)
        inp[:] = input - mean[:,np.newaxis]
        C = np.dot(inp,np.transpose(inp))
        self.s, self.v = np.linalg.eig(C)


    def project(self, inp, n=0):
        print('Projecting onto PCs...')
        output = np.dot(np.transpose(self.v),inp)
        return output

    def whiten(self,inp,eps=0.001):
        print('Whitening data...')
        # whiten the data
        d = np.sqrt(self.s)
        d = np.where(np.abs(d) > eps,d,eps)
        d = np.diag(1/np.sqrt(d))
        # print(d)
        mat = np.dot(self.v.T,inp)
        # mat = np.dot(self.v,mat)
        return np.dot(d,mat)

class ICA(PCA):
    """docstring for ICA."""

    def __init__(self, input, eps=0.001):
        pca = PCA(input)
        self.inp = pca.whiten(input,eps=eps)
        self.g = self.deriv_1
        self.gd = self.deriv_2


    def iterate(self,maxiter=10,eps=1e-12,ncomp=0):

        if ncomp == 0:
            nc = self.inp.shape[0]
        else:
            nc = ncomp

        self.w = np.zeros((nc,self.inp.shape[0]))
        for c in range(nc):
            w = np.random.rand(self.inp.shape[0])
            error = 1
            iter = 0
            while iter < maxiter and error > eps:
                old_w = w
                wTx = np.dot(old_w.T,self.inp)
                mat = self.g(wTx)
                mat = np.dot(self.inp,mat.T) / self.inp.shape[1]
                w = mat - old_w * np.mean(self.gd(wTx))

                # project orthogonal
                w = w - np.dot(np.dot(self.w[0:c,:],w).T,self.w[0:c,:])
                w = w / np.linalg.norm(w)

                error = np.abs(np.dot(w,old_w))-1
                iter = iter + 1
                print('N%d: iteration %d, error = %f' % (c,iter,error))

            self.w[c,:] = w

    def project(self,input):
        return np.dot(self.w,input)

    def fun(self,x):
        return -np.exp(-x**2/2)

    def deriv_1 (self,x):
        # return np.tanh(x)
        return x * np.exp(-x**2/2)

    def deriv_2(self,x):
        # return 1 - np.tanh(x)*np.tanh(x)
        return (1-x**2) * np.exp(-x**2/2)
