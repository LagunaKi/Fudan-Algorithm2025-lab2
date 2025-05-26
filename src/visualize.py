import argparse
import ast
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"

def read_output(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    final_ans = ast.literal_eval(lines[0].split('=',1)[1].strip())
    anchors = ast.literal_eval(lines[1].split('=',1)[1].strip())
    final_ans_strand = ast.literal_eval(lines[2].split('=',1)[1].strip())
    anchors_strand = ast.literal_eval(lines[3].split('=',1)[1].strip())
    return final_ans, anchors, final_ans_strand, anchors_strand

def visualize_alignment_plotly(final_ans, anchors, final_ans_strand, anchors_strand, query_len, ref_len):
    fig = go.Figure()
    # final_ans主链
    for (q0, q1, r0, r1), strand in zip(final_ans, final_ans_strand):
        color = 'blue' if strand == 1 else 'violet'
        name = '主链-正向' if strand == 1 else '主链-反向'
        fig.add_trace(go.Scatter(
            x=[q0, q1], y=[r0, r1],
            mode='lines',
            line=dict(color=color, width=3),
            name=name,
            legendgroup=name,
            showlegend=True
        ))
    # anchors（正向红色，反向半透明橙色，默认legendonly）
    if anchors:
        # 正向
        x_anchors_pos, y_anchors_pos = [], []
        x_anchors_neg, y_anchors_neg = [], []
        for (q0, q1, r0, r1), strand in zip(anchors, anchors_strand):
            if strand == 1:
                x_anchors_pos += [q0, q1, None]
                y_anchors_pos += [r0, r1, None]
            else:
                x_anchors_neg += [q0, q1, None]
                y_anchors_neg += [r0, r1, None]
        if x_anchors_pos:
            fig.add_trace(go.Scatter(
                x=x_anchors_pos, y=y_anchors_pos,
                mode='lines',
                line=dict(color='red', width=1, dash='dot'),
                name='锚点-正向',
                legendgroup='锚点-正向',
                opacity=1.0,
                showlegend=True,
                visible='legendonly'
            ))
        if x_anchors_neg:
            fig.add_trace(go.Scatter(
                x=x_anchors_neg, y=y_anchors_neg,
                mode='lines',
                line=dict(color='orange', width=1, dash='dot'),
                name='锚点-反向',
                legendgroup='锚点-反向',
                opacity=0.4,
                showlegend=True,
                visible='legendonly'
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
        legend_title='图例',
    )
    fig.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA比对结果可视化工具（plotly版）')
    parser.add_argument('--input', type=str, default='output.txt', help='输入output.txt')
    args = parser.parse_args()

    final_ans, anchors, final_ans_strand, anchors_strand = read_output(args.input)
    # 自动推断query/ref长度
    all_quads = final_ans + anchors
    query_len = max([max(q0, q1) for q0, q1, _, _ in all_quads]) + 1 if all_quads else 0
    ref_len = max([max(r0, r1) for _, _, r0, r1 in all_quads]) + 1 if all_quads else 0
    visualize_alignment_plotly(final_ans, anchors, final_ans_strand, anchors_strand, query_len, ref_len)

    '''
    python src/visualize.py --input output.txt
    python src/visualize.py --input output.txt --img_size 3000 --query_len 10000 --ref_len 8000
    '''