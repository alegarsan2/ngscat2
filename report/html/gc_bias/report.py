import numpy as np
from plotly import graph_objs as go
import plotly
import scipy.stats as st


class Report():
    def __init__(self,outdir):
        self.outdir = outdir
    def report(self, gclist, meanlists, coveragefiles):
        ymax = max(max(meanlist) for meanlist in meanlists)
        ymin = min(min(meanlist) for meanlist in meanlists)
        xmax = max(gclist)
        xmin = min(gclist)
        colorsc = [[0.0, '#fcf9f7'],
                    [0.16666666666666666, '#edcfc9'],
                    [0.3333333333333333, '#daa2ac'],
                    [0.5, '#bc7897'],
                    [0.6666666666666666, '#925684'],
                    [0.8333333333333333, '#5f3868'],
                    [1.0, '#2d1e3e']]

        for indx, meanlist in enumerate(meanlists):
            self.make_kdeplot(gclist,meanlist, xmin, xmax, ymin, ymax, 100,
                              colorsc, title= coveragefiles[indx].name.decode('utf-8').split("/")[-1])





    def kde_scipy(self, vals1, vals2, a, b, c, d, N):
        # vals1, vals2 are the values of two variables (columns)
        # (a,b) interval for vals1; usually larger than (np.min(vals1), np.max(vals1))
        # (c,d) -"-          vals2

        x = np.linspace(a, b, N)
        y = np.linspace(c, d, N)
        X, Y = np.meshgrid(x, y)
        positions = np.vstack([Y.ravel(), X.ravel()])

        values = np.vstack([vals1, vals2])
        kernel = st.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)

        return [x, y, Z]

    def make_kdeplot(self, varX, varY, a, b, c, d, N, colorsc, title):
        # varX, varY are lists, 1d numpy.array(s), or dataframe columns, storing the values of two variables

        x, y, Z = self.kde_scipy(varY, varX, a, b, c, d, N)

        trace = [go.Contour(
                z=Z,
                x=x,
                y=y,
                colorscale = colorsc,
                # reversescale=True,
                opacity=0.9,
                contours= dict(showlines = False)
            )]

        # trace = [go.Mesh3d(
        #         z=Z,
        #         x=x,
        #         y=y,
        #         #colorscale = colorsc,
        #         # reversescale=True,
        #         opacity=0.9,
        #         #contours= dict(showlines = False)
        #     )]

        layout_comp = go.Layout(
            title=title,
            showlegend=False,
            autosize=True,
            # width=650,
            # height=650,
            xaxis=dict(
                title = 'GC content (%)',
                range=[a, b],
                showgrid=False,
                nticks=7
            ),
            yaxis=dict(
                title= 'Mean Coverage',
                range=[c, d],
                showgrid=False,
                nticks=7
            ),
            margin=go.layout.Margin(
                l=50,
                r=40,
                b=85,
                t=100,
            ),
        )

        fig = go.Figure(data=trace, layout=layout_comp)
        plotly.offline.plot(fig, filename=self.outdir + 'gcbias_plot.html',
                            auto_open=True,
                            config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))