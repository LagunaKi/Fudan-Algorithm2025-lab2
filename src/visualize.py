import argparse
import ast
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"

def read_output(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    # 解析两行
    final_ans = ast.literal_eval(lines[0].split('=',1)[1].strip())
    anchors = ast.literal_eval(lines[1].split('=',1)[1].strip())
    return final_ans, anchors

def visualize_alignment_plotly(final_ans, anchors, query_len, ref_len):
    fig = go.Figure()
    # final_ans主链（蓝色）
    for q0, q1, r0, r1 in final_ans:
        fig.add_trace(go.Scatter(
            x=[q0, q1], y=[r0, r1],
            mode='lines',
            line=dict(color='blue', width=2),
            name='final_ans',
            showlegend=True,
            visible=True
        ))
    # anchors（红色，所有anchors一次性绘制，默认隐藏）
    if anchors:
        x_anchors = []
        y_anchors = []
        for q0, q1, r0, r1 in anchors:
            x_anchors += [q0, q1, None]
            y_anchors += [r0, r1, None]
        fig.add_trace(go.Scatter(
            x=x_anchors, y=y_anchors,
            mode='lines',
            line=dict(color='red', width=1, dash='dot'),
            name='anchors',
            showlegend=True,
            visible='legendonly'  # 默认不显示，点击图例可一键显示/隐藏
        ))
    fig.update_layout(
        title='DNA序列比对可视化',
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
    parser = argparse.ArgumentParser(description='DNA比对结果可视化工具（plotly版）')
    parser.add_argument('--input', type=str, default='output.txt', help='输入output.txt')
    args = parser.parse_args()

    final_ans, anchors = read_output(args.input)
    # 自动推断query/ref长度
    query_len = max([max(q0, q1) for q0, q1, _, _ in final_ans+anchors]) + 1
    ref_len = max([max(r0, r1) for _, _, r0, r1 in final_ans+anchors]) + 1
    visualize_alignment_plotly(final_ans, anchors, query_len, ref_len)

    '''
    python src/visualize.py --input output.txt
    python src/visualize.py --input output.txt --img_size 3000 --query_len 10000 --ref_len 8000
    '''