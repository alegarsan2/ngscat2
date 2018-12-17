import plotly.graph_objs as go
import plotly
class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, coveragefiles, results):
        self.target_coverage_plot(coveragefiles, results)


    def target_coverage_plot(self, coveragefiles, target_coverage_results):
        colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
        data = []
        for i, result in enumerate(target_coverage_results):
            trace = go.Bar(
                x=list(result['perccoveredposition'].keys()),
                y=list(result['perccoveredposition'].values()),
                hoverinfo='text',
                hoverlabel=dict(font=dict(color=['black'])),
                text=['Percentage:' + str(value) + '%' for value in list(result['perccoveredposition'].values())],
                # mode='lines',
                name=coveragefiles[i].name.decode('utf-8').split('/t')[-1],
                marker=dict(color=colors[i],
                            line=dict(
                                color='rgb(0,0,0)',
                                width=.6)
                            ), )

            data.append(trace)

        layout_comp = go.Layout(
            title='% on positions covered',
            hovermode='closest',
            barmode='group',
            xaxis=dict(showticklabels=True, showgrid=True, title='Coverage threshold'),
            yaxis=dict(title='% covered positions', range=[0, 100]),
            margin=go.layout.Margin(
                l=50,
                r=10,
                b=10,
                t=50,
                pad=4
            ),
        )

        fig = go.Figure(data=data, layout=layout_comp)
        plotly.offline.plot(fig, filename=self.outdir + 'covered_positions.html',
                            auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))

