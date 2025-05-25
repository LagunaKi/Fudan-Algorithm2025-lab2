import argparse
import ast
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"
# 读取output.txt格式的匹配区间
def read_intervals(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        data = f.read().strip()
    # 解析为列表
    intervals = ast.literal_eval(data)
    return intervals

# plotly可视化主函数
def visualize_alignment_plotly(intervals, query_len, ref_len):
    fig = go.Figure()
    # 绘制每个匹配区间
    for q0, q1, r0, r1 in intervals:
        color = 'blue' if q0 <= q1 else 'red'
        fig.add_trace(go.Scatter(
            x=[q0, q1],
            y=[r0, r1],
            mode='lines',
            line=dict(color=color, width=2),
            showlegend=False
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

    intervals = read_intervals(args.input)
    # 自动推断query/ref长度
    query_len = max(max(q0, q1) for q0, q1, _, _ in intervals) + 1
    ref_len = max(max(r0, r1) for _, _, r0, r1 in intervals) + 1
    visualize_alignment_plotly(intervals, query_len, ref_len)

    '''
    python src/visualize.py --input output.txt
    python src/visualize.py --input output.txt --img_size 3000 --query_len 10000 --ref_len 8000
    '''