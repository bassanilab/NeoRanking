import sys
import numpy as np
from scipy.optimize import basinhopping
from matplotlib import pyplot as plt

class CalcUmapClusterProbThresh:

    def __init__(self, rank_file):

        self.rank_file = rank_file

        f = open(rank_file, 'r')
        lines = f.readlines()
        f.close()

        self.data = []
        self.parse(lines)

        rho_min = 0.0
        rho_max = 0.01
        alpha_min = - 1.0
        alpha_max = 1.0

        class Bounds(object):
            def __init__(self, xmax=[rho_max, alpha_max], xmin=[rho_min, alpha_min]):
                self.xmax = np.array(xmax)
                self.xmin = np.array(xmin)

            def __call__(self, **kwargs):
                x = kwargs["x_new"]
                tmax = bool(np.all(x <= self.xmax))
                tmin = bool(np.all(x >= self.xmin))
                return tmax and tmin

        x0 = [0.0,0.0]
        minimizer_kwargs = {"method": "BFGS"}
        print("initial value: x = [rho, alpha] = %s, rank_scoee(x) = %.4f" % (str(x0), self.calc_score(x0)))
        ret = basinhopping(self.calc_score, x0, minimizer_kwargs=minimizer_kwargs, niter=500, accept_test=Bounds())
        print("global minimum: x = [rho, alpha] = %s, rank_score(x) = %.4f" % (str(ret.x), ret.fun))

        fig = plt.figure(figsize=(10, 10))

        xx, yy = np.mgrid[rho_min:rho_max:40j, alpha_min:alpha_max:40j]
        ax = fig.gca()
        ax.set_xlim(rho_min, rho_max)
        ax.set_ylim(alpha_min, alpha_max)

        positions = np.vstack([xx.ravel(), yy.ravel()])

        f = np.zeros(positions.shape[1])
        for i in range(positions.shape[1]):
            f[i] = self.calc_score(positions[:,i])

        f = np.reshape(f, xx.shape)
#        cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
        ax.imshow(np.rot90(f), cmap='coolwarm', extent=[rho_min, rho_max, alpha_min, alpha_max], aspect='auto')
        cset = ax.contour(xx, yy, f, colors='k',levels=20)
        ax.plot(ret.x[0], ret.x[1], color='red', marker='*', linewidth = 2, markersize = 12)
        ax.clabel(cset, inline=1, fontsize=10)
        ax.set_xlabel('rho')
        ax.set_ylabel('alpha')
        plt.title('Best match [rho, alpha] = {0}, immuno_missed_top20(x) = {1:.0f}'.format(str(ret.x), ret.fun))

        # x = np.arange(-5, 10, 0.02).tolist()
        # y = np.zeros(len(x))
        # for i in range(len(x)):
        #     y[i] = self.calc_score([ret.x[0],x[i]])
        #
        #
        # plt.scatter(x,y)
        # plt.xlabel("alpha")
        # plt.ylabel("Fitness")

        plt.show()



    def parse(self, lines):

        count = 0
        # Strips the newline character
        m = np.ndarray(shape=[0,0])
        for line in lines:
            if count % 6 == 0:
                if count > 0:
                    self.data.append(m)
                idx = line.split(',')
                m = np.ndarray(shape=[5, len(idx)])
            else:
                m[(count-1)%6,:] = line.split(',')

            count += 1

    def calc_score(self, x):

        score = 0.0
        for m in self.data:
            y = np.array(m[0,:])
            y_pred = np.array(m[1,:])
            p_umap = np.array(m[4,:])

            y_adj = np.array(list(map(lambda y, p: 0 if p < x[0] else y * np.power(p, x[1]), y_pred, p_umap)))
#            y_adj = np.array(list(map(lambda y, p: y/2 if p<x[0] else y*np.power(1+p,x[1]), y_pred, p_umap)))
#            y_adj = np.array(list(map(lambda y, p: y*np.power(1+p,x[1]), y_pred, p_umap)))

            score += self.top_20(y,y_adj)

        return 56-score


    def sum_rank_correct(self, y_true, y_pred):
        idx = np.argsort(-y_pred)
        y_true = y_true[idx]

        r = np.where(y_true == 1)

        return np.sum(1-np.exp(np.multiply(-0.05, r)))


    def top_20(self, y_true, y_pred):
        idx = np.argsort(-y_pred)
        y_true = y_true[idx]

        return np.sum(y_true[:20]==1)

if __name__ == '__main__':

    CalcUmapClusterProbThresh(sys.argv[1])
