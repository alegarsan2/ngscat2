import plotly.graph_objs as go
import plotly
import numpy as np
class Report():
    #FIXME como puedo acceder al outdir del  main report sin pasarle el mainreport al report como atributo o como entrada.
    def __init__(self, mainreporter):
        self.mainreporter = mainreporter

        #Aquí añadiremos las keys necesarias para generar el report
        self.summary = {}
        self.plot_dir_hist= mainreporter.outdir + '/data/coverage_histplot.html'
        self.plot_dir_boxplot= mainreporter.outdir + '/data/coverage_boxplot.html'

    def report(self, coveragefiles, target_distribution_results):
        self.target_distribution_histplot(coveragefiles, target_distribution_results)
        self.target_distribution_boxplot(coveragefiles, target_distribution_results)

        warnthreshold = []
        status = []
        mean = []

        for i, result in enumerate(target_distribution_results):
            warnthreshold.append(result['warnthreshold'])
            status.append(result['status'])
            mean.append(result['mean'])

        self.summary['warnthreshold'] = warnthreshold
        self.summary['status'] = status
        self.summary['mean'] = mean

        self.mainreporter.addsection('distribution', self)



    def getsummary(self):
         #Attribute encapsulation
        return self.summary




    def target_distribution_histplot(self, coveragefiles, target_distribution_results):
        colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
        data = []
        for i, result in enumerate(target_distribution_results):
            trace = go.Bar(
                y=result['histdata']['numberread'],
                x=result['histdata']['coveragepos'],
                width=result['histdata']['width'],
                hoverinfo='text',
                hoverlabel=dict(font=dict(color=['black'])),
                text=['Count: ' + str(value) + ('<br>') + 'Coverage: ' + str(result['histdata']['coveragepos'][i])
                      for i, value in enumerate(result['histdata']['numberread'])],
                # mode='lines',
                opacity=0.7,
                name=coveragefiles[i].name.decode('utf-8').split('/')[-1],
                marker=dict(color=colors[i],
                            line=dict(
                                color='rgb(0,0,0)',
                                width=.6)
                            ),)
            data.append(trace)

        layout_comp = go.Layout(
            title='Histogram',
            hovermode='closest',
            barmode='group',
            xaxis=dict(showticklabels=True, showgrid=True, title='Coverage'),
            yaxis=dict(title='Count'),
            margin=go.layout.Margin(
                l=50,
                r=10,
                b=20,
                t=50,
                pad=8),
            legend= dict(bgcolor= 'rgba(0,0,0,0)')
        )

        fig = go.Figure(data=data, layout=layout_comp)
        plotly.offline.plot(fig, filename= self.plot_dir_hist,
                            auto_open=False, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))

    def target_distribution_boxplot(self, coveragelist, target_distribution_result):
        colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
        data = []
        for i, coveragefile in enumerate(coveragelist):
            trace = go.Box(
                # Random subsampling in order to represent fasther the data. Final size 100000 points
                y= np.random.choice(coveragefile.coverages,
                size=int(len(coveragefile.coverages) / ((len(coveragefile.coverages) // 100000)
                            if len(coveragefile.coverages) > 100000 else 1))),
                name= coveragefile.name.decode('utf-8').split('/')[-1],
                marker=dict(
                    color=colors[i],
                ),
                boxpoints='suspectedoutliers',
                jitter=0.01)
            data.append(trace)

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
                b=20,
                t=50,
                pad=8
            ),
        )
        fig = go.Figure(data=data, layout=layout_comp)
        plotly.offline.plot(fig, filename=self.plot_dir_boxplot,
                            auto_open=False, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'],
                                                         showlink=False))