import plotly.graph_objs as go
import plotly
import numpy as np

class Report():
    def __init__(self, mainreporter):

        self.plot_dir_hist = mainreporter.outdir + '/data/std_histplot.html'
        self.plot_dir_boxplot = mainreporter.outdir + '/data/std_boxplot.html'
        self.summary = {}
        self.mainreporter = mainreporter

    def report(self, stdlists, results):
        self.std_distr_box(stdlists, results)
        self.std_distr_hist(stdlists, results)


        status = []
        mean = []
        warnthreshold = []

        for i, result in enumerate(results['results']):
            # item = dict()
            # item['status'] = result['status']
            # item['mean'] = result['mean']
            #
            # self.summary[result['bamfilename']] = item

            status.append(result['status'])
            mean.append(result['mean'])


        self.summary['status']= status
        self.summary['mean'] = mean
        self.summary['warnthreshold'] = results['warnthreshold']
        self.mainreporter.addsection('stdev', self)

    def getsummary(self):
        # Attribute encapsulation
        return self.summary

    def std_distr_hist(self, stdlists, results):
        # Plot His
        traces = []

        for indx, stdlist in enumerate(stdlists):
            trace = go.Histogram(
                x=stdlist,
                xbins=dict(end=1, size=0.007, start=0),
                # hoverinfo='text',
                # hoverlabel=dict(font=dict(color=['black'])),
                # text=['Count: ' + str(x[i]) + '<br>' + 'Coverage: ' + str(results['histdata']['stdpos'][i])
                # for i, value in enumerate(results['histdata']['numberstd'])],
                # mode='lines',
                opacity=0.7,
                name=results['results'][indx]['bamfilename'].split('/')[-1],
                marker=dict(
                    line=dict(
                        color='rgb(0,0,0)',
                        width=.6)))
            traces.append(trace)

        layout_comp = go.Layout(
            title='Histogram',
            # hovermode='closest',
            barmode='overlay',
            xaxis=dict(showticklabels=True, showgrid=True, title='Normalized standart deviation', zeroline=True),
            yaxis=dict(showticklabels=True, title='Frequency', zeroline=True),
            showlegend=True,
            legend=dict(x=0.75,
                        y=1),
            margin=go.layout.Margin(
                l=60,
                r=20,
                b=30,
                t=50,
                pad=4
            ),
        )

        fig = go.Figure(data=traces, layout=layout_comp)
        plotly.offline.plot(fig, filename=self.plot_dir_hist,
                            auto_open=False,  show_link= False, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'],
                                                         showlink=False))


    def std_distr_box(self,stdlists, results):
        # Plot
        traces =[]
        for indx, stdlist in enumerate(stdlists):
            trace = go.Box(
                # Random subsampling in order to represent fasther the data. Final size 100000 points
                y=np.random.choice(stdlist, size=int(len(stdlist) / ((len(stdlist) // 100000) if
                                                                              len(stdlist) > 100000 else 1))),
                name=results['results'][indx]['bamfilename'].split('/')[-1],
                boxpoints='suspectedoutliers',
                jitter=0.01)
            traces.append(trace)

        layout_comp = go.Layout(
            title='',
            hovermode='closest',
            # barmode='group',
            xaxis=dict(showticklabels=True, showgrid=True, title=''),
            yaxis=dict(title='Depth',
                       autorange=True),
            margin=go.layout.Margin(
                l=50,
                r=10,
                b=10,
                t=50,
                pad=4
            ),
        )

        fig = go.Figure(data=traces, layout=layout_comp)
        plotly.offline.plot(fig, filename= self.plot_dir_boxplot,
                            auto_open=False, show_link= False,
                            config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))

#TODO esta parte es la que va en el main para llamar a a esta funcion.
# from 'report/html/depth_stdev/' import Reporter
# from depth_stdev_computer
#
# reporter = Report('/outdir')
# depth_stdev_computer(coveragefiles, reporter.report)
