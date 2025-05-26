'''
python std_2.py
'''
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"

# 标准区间列表
final_ans = [
    (0, 300, 0, 300),
    (300, 400, 400, 500),
    (400, 501, 499, 600),
    (501, 801, 601, 901),
    (801, 900, 701, 800),
    (900, 1000, 700, 800),
    (1000, 1200, 700, 900),
    (1200, 1300, 900, 1000),
    (1300, 1400, 900, 1000),
    (1402, 1501, 402, 501),
    (1501, 1618, 1001, 1118),
    (1618, 1700, 1318, 1400),
    (1700, 1802, 1200, 1302),
    (1807, 1900, 1107, 1200),
    (1900, 1998, 1400, 1498),
    (2299, 2500, 1499, 1700)
]

def visualize_alignment_plotly(final_ans):
    fig = go.Figure()
    # final_ans主链（蓝色）
    for q0, q1, r0, r1 in final_ans:
        fig.add_trace(go.Scatter(
            x=[q0, q1], y=[r0, r1],
            mode='lines',
            line=dict(color='green', width=2),
            name='final_ans',
            showlegend=True,
            visible=True
        ))
    fig.update_layout(
        title='标准DNA序列比对区间可视化',
        xaxis_title='Query',
        yaxis_title='Reference',
        width=800,
        height=800,
        plot_bgcolor='white',
        xaxis=dict(showgrid=True, gridcolor='lightgray', zeroline=False),
        yaxis=dict(showgrid=True, gridcolor='lightgray', zeroline=False),
    )
    fig.show()

if __name__ == '__main__':
    visualize_alignment_plotly(final_ans) 