import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from matplotlib import cm 
import networkx as nx
np.random.seed(42)

# # Page configure
# st.set_page_config(page_title="TanLab", layout="wide")

# # ===========================================================================================
# # 生成模拟数据
# drugs = ['DrugA', 'DrugB', 'DrugC', 'DrugD', 'DrugE']
# data = []

# for drug in drugs:
#     for i in range(1, 101):
#         data.append([drug, f'Target{i}', np.random.uniform(0, 1), np.random.uniform(-2, 2), np.random.uniform(0, 0.07)])

# df = pd.DataFrame(data, columns=['Drugs', 'Targets', 'Hscore', 'logFC', 'Pvalue'])
# ===========================================================================================

df = pd.read_excel("./data/Top100.xlsx")
print(df.columns)
col_names = ['Targets', 'Gene', 'logFC', 'logPvalue', 'logAdjPvalue', 'Class', 'Hscore', 'Drugs']
df.columns = col_names
df['UniprotURL'] = [f'https://www.uniprot.org/uniprotkb/{item}/entry' for item in df['Targets']]

print(df.head(5))

# 假设你有一个字典，存储文件名和文件路径
files = {
    "primary-screen-readme.txt": "./Data/primary-screen-readme.txt",
    "primary-screen-cell-line-info.csv": "./Data/primary-screen-cell-line-info.csv",
    "primary-screen-pooling-info.csv": "./Data/primary-screen-pooling-info.csv"}


# 顶部导航栏
tab1, tab2 = st.tabs(["🏠 Home", "🌎 Contact"])
with tab1:
    # 使用列布局来设置文章标题、作者和实验室logo
    with st.container():
        col1, col2, col3 = st.columns([1, 6, 1])
        with col1:
            pass
        #     tab1, tab2 = st.tabs(["Home", "Concat"])
            # st.image("https://faculty.sustech.edu.cn/wp-content/uploads/2019/11/2019110717454716.jpg", width=200)


        with col2:
            # 中间列放置标题和作者
            st.markdown(
                "<h1 style='text-align: center; color: black; background-color: white;'>Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling</h1>",
                unsafe_allow_html=True)

            st.markdown("""
            <p style='text-align: center; color: black; background-color: white;'>
            Corsello SM, Nagari RT, Spangler RD, Rossen J, Kocak M, Bryan JG, Humeidi R, Peck D, Wu X, Tang AA, Wang VM, Bender SA, Lemire E, Narayan R, Montgomery P, Ben-David U, Garvie CW, Chen Y, Rees MG, Lyons NJ, McFarland JM, Wong BT, Wang L, Dumont N, O'Hearn PJ, Stefan E, Doench JG, Harrington CN, Greulich H, Meyerson M, Vazquez F, Subramanian A, Roth JA, Bittker JA, Boehm JS, Mader CC, Tsherniak A, Golub TR
            </p>
            """, unsafe_allow_html=True)

        with col3:
            # 右侧放置实验室logo
            # image = Image("C:/Users/17608/Desktop/MyAPP/figures/logo.jpg")
            st.image("images/logo.jpg", width=200)

    # 摘要部分
    st.subheader("Abstract")
    st.write("""
    Anti-cancer uses of non-oncology drugs have occasionally been found, but such discoveries have been serendipitous. We sought to create a public resource containing the growth inhibitory activity of 4,518 drugs tested across 578 human cancer cell lines. We used PRISM, a molecular barcoding method, to screen drugs against cell lines in pools. An unexpectedly large number of non-oncology drugs selectively inhibited subsets of cancer cell lines in a manner predictable from the cell lines' molecular features. Our findings include compounds that killed by inducing PDE3A-SLFN12 complex formation; vanadium-containing compounds whose killing depended on the sulfate transporter SLC26A2; the alcohol dependence drug disulfiram, which killed cells with low expression of metallothioneins; and the anti-inflammatory drug tepoxalin, which killed via the multi-drug resistance protein ABCB1. The PRISM drug repurposing resource ([link](https://depmap.org/repurposing)) is a starting point to develop new oncology therapeutics, and more rarely, for potential direct clinical translation.
    """)

    landscape, network = st.tabs(["1️⃣ Scatter Plot", "2️⃣ Network Interaction"])

    with landscape:
        col1, col2, col3 = st.columns([1, 0.1, 1])
        selected_drug = col1.selectbox(":four_leaf_clover: Select a Drug", df['Drugs'].unique())
        selected_hscore = col3.slider(':herb: Set Hscore threshold', min_value=0.0, max_value=1.0, value=0.8, step=0.01)
        filtered_df = df[(df['Drugs'] == selected_drug) & (df['Hscore'] > selected_hscore)]

        # 绘制火山图
        col1, col2, col3 = st.columns([2, 5, 2])
        if not filtered_df.empty:
            fig = px.scatter(filtered_df, 
                                x='logFC', 
                                y='logPvalue', 
                                color='Hscore', 
                                labels={'x': 'logFC', 'y': 'logPvalue', 'target':'Targets'},
                                # color_continuous_scale=px.colors.diverging.RdBu,
                                color_continuous_scale="reds",
                                hover_data={'Targets': True}
                                )
            fig.update_traces(marker_size=12)
            fig.update_layout(
            margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
            autosize=False,
            width=1000,   # 设置图表的宽度
            height=500,
            title_x=0.5,  # 设置标题居中
            xaxis_title_font=dict(size=20, weight='bold'),  # x轴标题字体大小
            yaxis_title_font=dict(size=20, weight='bold'),  # y轴标题字体大小
            xaxis_tickfont=dict(size=18),  # x轴刻度字体大小
            yaxis_tickfont=dict(size=18),)
            # xaxis=dict(showgrid=True, zeroline=True, showticklabels=True),
            # yaxis=dict(showgrid=True, zeroline=True, showticklabels=True))  # y轴刻度字体大小
        
            # fig.add_shape(type="line",
            #         x0=0, y0=0, x1=0, y1=max(filtered_df['logPvalue']),
            #         line=dict(color='#EAEAEA', width=2, dash="dot"),
            #         name="x=0 line")
            
            col2.plotly_chart(fig, use_container_width=True)
        else:
            col2.write('No results found.')

    # with biomarker:
    #     # reference ----> https://plotly.com/python/box-plots/
    #     swarm_df = px.data.tips()
    #     fig = px.strip(swarm_df, x='day', y='tip')
    #     fig.update_layout(
    #         margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
    #         autosize=False,
    #         width=500,  # 你可以设置图表的宽度
    #         height=400,
    #         title_x=0.5,  # 设置标题居中
    #         xaxis_title_font=dict(size=20),  # x轴标题字体大小
    #         yaxis_title_font=dict(size=20),  # y轴标题字体大小
    #         xaxis_tickfont=dict(size=18),  # x轴刻度字体大小
    #         yaxis_tickfont=dict(size=18)  # y轴刻度字体大小
    #     )
    #     st.plotly_chart(fig, use_container_width=True)

    with network:
        col11, col22, col33 = st.columns([1, 1, 1])
        selected_drug = col11.selectbox(":maple_leaf: Select a Drug", df['Drugs'].unique())
        layout= col22.selectbox(':fallen_leaf: Choose a network layout',('Random Layout','Spring Layout','Shell Layout','Kamada Kawai Layout'))
        # selected_hscore = col33.text_input('Set Hscore threshold')
        selected_hscore = col33.slider(':leaves: Set Hscore threshold', min_value=0.0, max_value=1.0, value=0.8, step=0.01)
        filtered_df = df[(df['Drugs'] == selected_drug) & (df['Hscore'] > selected_hscore)]

        # 绘制火山图
        col1, col2, col3 = st.columns([0.5, 5, 1])
        if not filtered_df.empty:
            
            # Create graph
            G = nx.Graph()

            # 添加节点和边 'Drugs', 'Target', 'Hscore', 'logFC', 'logPvalue'
            
            for index, row in filtered_df.iterrows():
                if row['Drugs'] not in G.nodes():
                    G.add_node(row['Drugs'], type='drug')  # Drug节点添加颜色和大小
                if row['Targets'] not in G.nodes():
                    G.add_node(row['Targets'], type='target', score=row['Hscore'])  # Target节点为圆形
                if not G.has_edge(row['Drugs'], row['Targets']):
                    G.add_weighted_edges_from([(row['Drugs'], row['Targets'], row['Hscore'])])
                    
                    # G.add_edge(row['Drugs'], row['Targets'], weight=row['Hscore'])  
                    # G.add_edge(row['Drugs'], row['Targets'], weight=row['Hscore'])
            
            # 获取节点位置
            #Get the position of each node depending on the user' choice of layout
            if layout=='Random Layout':
                pos = nx.random_layout(G) 
            elif layout=='Spring Layout':
                pos = nx.spring_layout(G, k=0.5, iterations=50)
            elif layout=='Shell Layout':
                pos = nx.shell_layout(G)            
            elif layout=='Kamada Kawai Layout':
                pos = nx.kamada_kawai_layout(G) 
            
            #Add positions of nodes to the graph
            for n, p in pos.items():
                G.nodes[n]['pos'] = p

            # 创建边的trace
            edge_trace = go.Scatter(x=[], 
                                    y=[], 
                                    line=dict(width=2, color='#EFEFEF'), ## '#EFEFEF'
                                    hoverinfo='none', 
                                    mode='lines')
            
            weights = [data['weight'] for _, _, data in G.edges(data=True)]
            min_weight = min(weights)
            max_weight = max(weights)
            cmap = cm.get_cmap('Greens')


            for edge in G.edges(data=True):
                x0, y0 = G.nodes[edge[0]]['pos']
                x1, y1 = G.nodes[edge[1]]['pos']
                edge_trace['x'] += tuple([x0, x1, None])
                edge_trace['y'] += tuple([y0, y1, None])
            


            # Adding nodes to plotly scatter plot
            node_trace = go.Scatter(x=[], 
                                    y=[], 
                                    text=[], 
                                    mode='markers+text', 
                                    textposition="top center", 
                                    hoverinfo='text',
                                    marker=dict(
                                                showscale=True,
                                                size=[],
                                                color=[],
                                                symbol =[],
                                                colorscale='blues',  # RdBu
                                                colorbar=dict(
                                                            title="Hscore",  # 颜色条的标题
                                                            titleside="top",  # 标题在右边
                                                            # tickmode="array",  # 自定义刻度
                                                            # tickvals=[0, 1],  # 刻度值，例如最低和最高的Hscore
                                                            # ticktext=["Low", "High"],  # 对应的标签
                                                            # ticks="outside"  # 刻度显示在外侧
                                                        )
                                                ))
                        
            
            for node in G.nodes():
                x, y = G.nodes[node]['pos']
                node_trace['x'] += tuple([x])
                node_trace['y'] += tuple([y])
                node_trace['text'] += tuple([node])  # Node labels
                node_type = G.nodes[node]['type']
                if node_type == 'drug':
                    node_trace['marker']['symbol'] += tuple(['star'])   # Star shape for drugs
                    node_trace['marker']['color'] += tuple(['red'])  # Red color for drug nodes
                    node_trace['marker']['size'] += tuple([25])  # Fixed size for drug nodes
                elif node_type == 'target':
                    # node_trace['marker']['color']  += tuple(['#7E2678'])   # 靶点节点蓝色#53AAA9
                    # node_trace['marker']['symbol'] = 'circle'  # Circle shape for targets
                    node_trace['marker']['color'] += tuple([G.nodes[node]['score']])  # Color based on Hscore
                    node_trace['marker']['size'] += tuple([25])  # Fixed size for target nodes
                


            for node, adjacencies in enumerate(G.adjacency()):
                # print('-------->', tuple([len(adjacencies[1])]))
                # node_trace['marker']['color'] += tuple([len(adjacencies[1])])  # Coloring each node based on the number of connections 
                node_info = adjacencies[0]
                # print('-------->', tuple([node_info]))
                node_trace['text'] += tuple([node_info])

            # 创建Plotly图表
            fig = go.Figure(data=[edge_trace, node_trace], 
                            layout=go.Layout(title='', #title takes input from the user
                                                title_x=0.45,
                                                titlefont=dict(size=25),
                                                showlegend=False,
                                                hovermode='closest',
                                            #  margin=dict(b=20,l=5,r=5,t=40),
                                                # margin=dict(b=5,l=5,r=5,t=5),
                                                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

            fig.update_layout(
            # margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
            autosize=False,
            width=500,  # 你可以设置图表的宽度
            height=500)
            col2.plotly_chart(fig, use_container_width=True)  

        else:
            col2.write(f"No data available for {selected_drug} with Hscore > {selected_hscore}.")


    # --------------------------------------------------------------------------------------------------------------------------------------------------------
    st.divider()
    left_col, media_col, right_col = st.columns([1, 0.1, 1])

    with left_col:
        st.subheader("PRISM Repurposing dataset")
        st.write("""
            The SIGMA Repurposing release contains small molecule viability datasets generated using the Broad Repurposing Library and the PRISM multiplexed cell-line viability assay.

            The primary PRISM Repurposing dataset contains the results of pooled-cell line chemical-perturbation viability screens for 4,518 compounds screened against 578 or 562 cell lines.

            The secondary PRISM Repurposing dataset contains the results of pooled-cell line chemical-perturbation viability screens for 1,448 compounds screened against 499 cell lines in a 8-step, 4-fold dilution, starting from 10μM.

            Data processing steps are described in the README file. Further descriptions of methods will be published in an upcoming publication.
        """)

    with media_col:
        pass

    with right_col:
        st.header("Datasets")
        st.subheader("Primary screen")
        # 遍历字典，生成下载按钮
        for file_name, file_path in files.items():
            # 读取文件内容
            with open(file_path, 'rb') as f:
                file_data = f.read()

            # 创建下载按钮
            st.download_button(
                label=file_name,  # 按钮标签
                data=file_data,  # 文件内容
                file_name=file_name,  # 下载时的文件名称
                mime='text/csv' if file_name.endswith('.csv') else 'text/plain'  # 文件类型
            )

    # --------------------------------------------------------------------------------------------------------------------------------------------------
with tab2:
    # st.subheader("Contact Us")
    # st.subheader("Contact Us")
    st.write("""🎈 Please feel free to contact christan@sustech.edu.cn with any issues, comments, or questions.""")

    st.text_area("**Message**")
    st.text_input("**Your Email**")
    if st.button("**Send**"):
        st.success("Message sent! We will get back to you soon.")


 # 网页底部版权信息
st.markdown(
    '<div style="display: flex; justify-content: center; align-items: center; '
    'height: 100px; background-color: #F2F2F2; color: black;">'
    '<p style="margin: 0;">Copyright 1995-2024 Kanehisa Laboratories</p>'
    '</div>',
    unsafe_allow_html=True
)