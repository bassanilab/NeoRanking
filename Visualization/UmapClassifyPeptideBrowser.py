import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import rankdata
from sklearn.cluster import DBSCAN
from UmapClusterProb import UmapClusterProb


class UmapClassifyPeptideBrowser:

    def __init__(self, umaper, df, y, y_pred, patients, peptides, name, plot_feature_directions = False):
        self.select = []

        self.umaper = umaper
        self.name = name
        self.y = y
        self.df = df
        self.x_umap = self.umaper.fill_missing_values(self.df.iloc[:, :])

        self.plot_directions = plot_feature_directions
        self.y_pred = y_pred
        self.patients = patients
        self.peptides = peptides

        self.x0_umap = 0
        self.x1_umap = 0
        self.y0_umap = 0
        self.y1_umap = 0
        self.artist = None
        self.N = self.df.shape[0]

        self.plot_umap()


    def plot_umap(self):

        self.plot_state = np.zeros(self.df.shape[0])

        self.plot_state[self.y == 1] = 2
        self.ranks = rankdata(-self.y_pred, method='min')

        for p in np.unique(self.patients):
            idx = np.where(self.patients == p)[0]
            idx_sort = np.argsort(-self.y_pred[idx])
            self.plot_state[idx[idx_sort[:20]]] = self.plot_state[idx[idx_sort[:20]]] + 1

        self.x_umap_min = np.min(self.x_umap[:, 0])
        self.x_umap_max = np.max(self.x_umap[:, 0])
        self.y_umap_min = np.min(self.x_umap[:, 1])
        self.y_umap_max = np.max(self.x_umap[:, 1])

        self.x_eps = abs(self.x_umap_max-self.x_umap_min)/100
        self.y_eps = abs(self.y_umap_max-self.y_umap_min)/100

        self.x0_umap = np.min(self.x_umap[:, 0])
        self.x1_umap = np.max(self.x_umap[:, 0])
        self.y0_umap = np.min(self.x_umap[:, 1])
        self.y1_umap = np.max(self.x_umap[:, 1])

        idx_neg_last = self.plot_state == 0
        idx_neg_first = self.plot_state == 1
        idx_cd8_last = self.plot_state == 2
        idx_cd8_first = self.plot_state == 3

        if self.name == "Rosenberg":
            plt.figure(figsize=(10, 10))
            plt.scatter(self.x_umap[:, 0], self.x_umap[:, 1], c=self.umaper.get_cluster_id(self.df.iloc[:,:]),
                        s=5, picker=True, label='cluster id', alpha=0.6, cmap='viridis')

        self.fig = plt.figure(figsize=(10, 10))

        self.ax1 = self.fig.add_subplot(2, 1, 1)

        if self.name != "Rosenberg":
            coords = self.umaper.fill_missing_values()
            self.ax1.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=5, picker=True, label='negative 20+',
                             alpha=0.3, cmap='viridis')

        self.ax1.scatter(self.x_umap[idx_neg_last, 0], self.x_umap[idx_neg_last, 1], c='lightgreen', s=5,
                         picker=True, label='negative 20+', alpha=0.5)
        self.ax1.scatter(self.x_umap[idx_neg_first, 0], self.x_umap[idx_neg_first, 1], c='green', s=5,
                         picker=True, label='negative top 20', alpha=0.7)
        self.ax1.scatter(self.x_umap[idx_cd8_first, 0], self.x_umap[idx_cd8_first, 1], c='red', s=5,
                         picker=True, label='CD8 top 20')
        self.ax1.scatter(self.x_umap[idx_cd8_last, 0], self.x_umap[idx_cd8_last, 1], c='orange', s=5,
                         picker=True, label='CD8 20+')

        plt.title("UMAP projection of {0} dataset".format(self.name), fontsize=10)
        plt.xlabel("UMAP dimension 1" )
        plt.ylabel("UMAP dimension 2" )

        lineHandle_1 = Line2D([], [], color='red', ls="-", label="CD8 & rank <= 20")
        lineHandle_2 = Line2D([], [], color='orange', ls="-", label="CD8 & rank > 20")
        lineHandle_3 = Line2D([], [], color='green', ls="-", label="Negative & rank <= 20")
        lineHandle_4 = Line2D([], [], color='lightgreen', ls="-", label="Negative & rank > 20")
        self.ax1.legend(handles=[lineHandle_1, lineHandle_2, lineHandle_3, lineHandle_4], loc=1, prop={'size': 5})

        self.ax2 = self.fig.add_subplot(2, 1, 2)

        self.fig.canvas.mpl_connect('button_press_event', lambda event: self.on_press_umap(event))
        self.fig.canvas.mpl_connect('button_release_event', lambda event: self.on_release_umap(event))
        self.fig.canvas.mpl_connect('key_press_event', lambda event: self.on_press(event))

        self.fig.show()


    def get_umaper(self):
        return self.umaper


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

        # RGBA (red, green, blue, alpha)
        col_1 = 'red'
        col_2 = 'orange'
        col_3 = [0.0, 0.3, 0.0, 0.7]
        col_4 = [0.0, 0.8, 0.0, 0.3]

        self.ax2.cla()

        if sum(self.select)>1:
            plt.title('Multiselection', fontsize=8)
        else:
            response = "CD8" if any(self.y[self.select]==1) else "negative"
            plt.title("{0},{1},{2:.3f},{3:.0f}".
                      format(str(np.array(self.peptides[self.select])[0]),
                             response,
                             np.array(self.y_pred[self.select])[0],
                             self.ranks[self.select][0]),
                      fontsize=8)

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


    def on_press_umap(self, event):
        self.x0_umap = event.xdata
        self.y0_umap = event.ydata


    def on_release_umap(self, event):

        self.x1_umap = event.xdata
        self.y1_umap = event.ydata

        if self.x0_umap > self.x1_umap:
            tmp = self.x0_umap
            self.x0_umap = self.x1_umap
            self.x1_umap = tmp

        if self.y0_umap > self.y1_umap:
            tmp = self.y0_umap
            self.y0_umap = self.y1_umap
            self.y1_umap = tmp

        self.select = (self.x_umap[:, 0] >= self.x0_umap-self.x_eps) & (self.x_umap[:, 1] >= self.y0_umap-self.y_eps) & \
                      (self.x_umap[:, 0] <= self.x1_umap+self.x_eps) & (self.x_umap[:, 1] <= self.y1_umap+self.y_eps)

        #        print(str(self.x0_umap) + ', ' + str(self.y0_umap) + ', ' + str(self.x1_umap) + ', ' + str(self.y1_umap))

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
