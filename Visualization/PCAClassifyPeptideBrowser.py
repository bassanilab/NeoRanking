import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from sklearn.decomposition import PCA
from scipy.stats import rankdata


class PCAClassifyPeptideBrowser:

    def __init__(self, df, y, y_pred, name, plot_feature_directions=False, interactive=False):
        self.select = []

        self.df = df
        self.name = name
        self.plot_directions = plot_feature_directions
        self.interactive = interactive
        self.y = y
        self.y_pred = y_pred

        self.x0_pca = 0
        self.x1_pca = 0
        self.y0_pca = 0
        self.y1_pca = 0
        self.x_pca = None
        self.artist = None
        self.N = df.shape[0]

        self.fig = plt.figure(figsize=(10, 10))

        if self.interactive:
            self.ax1 = self.fig.add_subplot(2, 1, 1)
        else:
            self.ax1 = self.fig.add_subplot(1, 1, 1)

        self.plot_state = np.zeros(df.shape[0])
        self.plot_state[y == 1] = 2
        self.ranks = rankdata(-y_pred, method='min')

        pca = PCA(n_components=2)

        self.x_pca = pca.fit_transform(df.iloc[:, :])
        variance = pca.explained_variance_ratio_

        self.x_pca_min = np.min(self.x_pca[:, 0])
        self.x_pca_max = np.max(self.x_pca[:, 0])
        self.y_pca_min = np.min(self.x_pca[:, 1])
        self.y_pca_max = np.max(self.x_pca[:, 1])

        self.x_eps = abs(self.x_pca_max-self.x_pca_min)/100
        self.y_eps = abs(self.y_pca_max-self.y_pca_min)/100

        self.x0_pca = np.min(self.x_pca[:, 0])
        self.x1_pca = np.max(self.x_pca[:, 0])
        self.y0_pca = np.min(self.x_pca[:, 1])
        self.y1_pca = np.max(self.x_pca[:, 1])

        idx_neg_last = self.plot_state == 0
        idx_neg_first = self.plot_state == 1
        idx_cd8_last = self.plot_state == 2
        idx_cd8_first = self.plot_state == 3

        self.ax1.scatter(self.x_pca[idx_neg_last, 0], self.x_pca[idx_neg_last, 1], c='lightgray', s=5,
                         picker=True, label='negative 20+', alpha=0.3)
        self.ax1.scatter(self.x_pca[idx_neg_first, 0], self.x_pca[idx_neg_first, 1], c='gray', s=5,
                         picker=True, label='negative top 20', alpha=0.7)
        self.ax1.scatter(self.x_pca[idx_cd8_first, 0], self.x_pca[idx_cd8_first, 1], c='red', s=5,
                         picker=True, label='CD8 top 20')
        self.ax1.scatter(self.x_pca[idx_cd8_last, 0], self.x_pca[idx_cd8_last, 1], c='blue', s=5,
                         picker=True, label='CD8 20+')

        plt.title("PCA projection of {0} dataset".format(self.name), fontsize=10)
        plt.xlabel("PC 1 (%.1f%%)" % (variance[0] * 100))
        plt.ylabel("PC 2 (%.1f%%)" % (variance[1] * 100))
        self.fig.tight_layout()

        legend = self.plot_feature_directions(pca)
        if legend is not None:
            plt.gca().add_artist(legend)

        lineHandle_1 = Line2D([], [], color='red', ls="-", label="CD8 & rank <= 20")
        lineHandle_2 = Line2D([], [], color='blue', ls="-", label="CD8 & rank > 20")
        lineHandle_3 = Line2D([], [], color='gray', ls="-", label="Negative & rank <= 20")
        lineHandle_4 = Line2D([], [], color='lightgray', ls="-", label="Negative & rank > 20")
        self.ax1.legend(handles=[lineHandle_1, lineHandle_2, lineHandle_3, lineHandle_4], loc=1, prop={'size': 5})

        if self.interactive:
            self.ax2 = self.fig.add_subplot(2, 1, 2)
        else:
            self.ax2 = None

        self.fig.canvas.mpl_connect('button_press_event', lambda event: self.on_press_pca(event))
        self.fig.canvas.mpl_connect('button_release_event', lambda event: self.on_release_pca(event))
        self.fig.canvas.mpl_connect('key_press_event', lambda event: self.on_press(event))

        self.fig.show()

    def plot_feature_directions(self, pca):

        if not self.plot_directions:
            return None

        coeff = np.transpose(pca.components_[0:2, :])
        n = coeff.shape[0]

        cmap = cm.get_cmap('hsv')
        handles = []
        for i in range(n):
            plt.arrow(0, 0, 3*coeff[i, 0], 3*coeff[i, 1], color=cmap((1.0*i)/n), width=0.01, head_width=0.06, alpha=0.5,
                      length_includes_head = True)
            handles.append(Line2D([], [], color=cmap((1.0*i)/n), ls="-", label=self.df.columns[i]))

        return self.ax1.legend(handles=handles, loc=4, prop={'size': 5})

    def plot_parallel_coords(self):

        if sum(self.select) == 0:
            return

        print("")
        print("**********************************************")
        print("Selection from {0}:".format(self.name))
        print("")
        print("Dataset\tPeptideId\tResponse\tScore\tRank\tTop20")
        n = self.df.shape[1]
        indices, = np.where(self.select)
        for i in indices:
            response = "CD8" if self.y[i]==1 else "negative"
            top20 = self.plot_state[i] in [1,3]
            print("{0}\t{1}\t{2}\t{3:.3f}\t{4:.0f}\t{5}".
                  format(self.name,str(self.patients[i]),response, self.y_pred[i], self.ranks[i],
                         str(top20)))

        df_1 = self.df.loc[np.logical_and(self.select, self.plot_state == 3),]
        df_1.insert(df_1.shape[1], "response", self.y[np.logical_and(self.select, self.plot_state == 3)])
        df_2 = self.df.loc[np.logical_and(self.select, self.plot_state == 2),]
        df_2.insert(df_2.shape[1], "response", self.y[np.logical_and(self.select, self.plot_state == 2)])
        df_3 = self.df.loc[np.logical_and(self.select, self.plot_state == 1),]
        df_3.insert(df_3.shape[1], "response", self.y[np.logical_and(self.select, self.plot_state == 1)])
        df_4 = self.df.loc[np.logical_and(self.select, self.plot_state == 0),]
        df_4.insert(df_4.shape[1], "response", self.y[np.logical_and(self.select, self.plot_state == 0)])

        col_1 = [1, 0, 0, 0.9]
        col_2 = [0, 0, 1, 0.7]
        col_3 = [0.3, 0.3, 0.3, 0.7]
        col_4 = [0.8, 0.8, 0.8, 0.3]

        self.ax2.cla()

        if df_4.shape[0] > 0:
            ax = pd.plotting.parallel_coordinates(df_4, class_column='response', cols=df_4.columns[0:n],
                                                  ax=self.ax2,
                                                  sort_labels=True,
                                                  color=[col_4])
            ax.xaxis.set_tick_params(labelsize=4, labelrotation=15)

        if df_3.shape[0] > 0:
            ax = pd.plotting.parallel_coordinates(df_3, class_column='response', cols=df_3.columns[0:n],
                                                  ax=self.ax2,
                                                  sort_labels=True,
                                                  color=[col_3])
            ax.xaxis.set_tick_params(labelsize=4, labelrotation=15)

        if df_2.shape[0] > 0:
            ax = pd.plotting.parallel_coordinates(df_2, class_column='response', cols=df_2.columns[0:n],
                                                  ax=self.ax2,
                                                  sort_labels=True,
                                                  color=[col_2])
            ax.xaxis.set_tick_params(labelsize=4, labelrotation=15)

        if df_1.shape[0] > 0:
            ax = pd.plotting.parallel_coordinates(df_1, class_column='response', cols=df_1.columns[0:n],
                                                  ax=self.ax2,
                                                  sort_labels=True,
                                                  color=[col_1])
            ax.xaxis.set_tick_params(labelsize=4, labelrotation=15)

        lineHandle_1 = Line2D([], [], color=col_1, ls="-", label="Immunogenic & rank <= 20")
        lineHandle_2 = Line2D([], [], color=col_2, ls="-", label="Immunogenic & rank > 20")
        lineHandle_3 = Line2D([], [], color=col_3, ls="-", label="Not immunogenic & rank <= 20")
        lineHandle_4 = Line2D([], [], color=col_4, ls="-", label="Not immunogenic & rank > 20")
        self.ax2.legend(handles=[lineHandle_1, lineHandle_2, lineHandle_3, lineHandle_4], loc=1, prop={'size': 6})

        plt.draw()


    def on_press_pca(self, event):
        self.x0_pca = event.xdata
        self.y0_pca = event.ydata


    def on_release_pca(self, event):

        self.x1_pca = event.xdata
        self.y1_pca = event.ydata

        if self.x0_pca > self.x1_pca:
            tmp = self.x0_pca
            self.x0_pca = self.x1_pca
            self.x1_pca = tmp

        if self.y0_pca > self.y1_pca:
            tmp = self.y0_pca
            self.y0_pca = self.y1_pca
            self.y1_pca = tmp

        self.select = (self.x_pca[:, 0] >= self.x0_pca-self.x_eps) & (self.x_pca[:, 1] >= self.y0_pca-self.y_eps) & \
                      (self.x_pca[:, 0] <= self.x1_pca+self.x_eps) & (self.x_pca[:, 1] <= self.y1_pca+self.y_eps)

#        print(str(self.x0_pca) + ', ' + str(self.y0_pca) + ', ' + str(self.x1_pca) + ', ' + str(self.y1_pca))

        self.plot_parallel_coords()


    def on_press(self, event):
        if len(self.select) == 0:
            return
        if event.key not in ('n', 'p'):
            return
        if event.key == 'n':
            inc = 1
        else:
            inc = -1

        self.select = np.array([x+inc if (0 <= x + inc < self.df.shape[0]) else -1 for x in self.select])
        self.select = self.select[self.select>=0]

        self.plot_parallel_coords()



