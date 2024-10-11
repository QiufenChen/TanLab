import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from matplotlib import cm 
import networkx as nx
np.random.seed(42)

# Page configure
st.set_page_config(page_title="📊", layout="wide")


my_df = pd.read_excel("./Data/OurData.xlsx")
idx = my_df.groupby(['Group', 'Drug', 'DrugName', 'Protein_ID', 'Gene'])['H-Score'].idxmax()
df_unique = my_df.loc[idx]
df = df_unique.groupby('Drug').filter(lambda x: len(x) >= 10)
df['fc'] = df.apply(lambda row: -abs(row['fc']) if row['lasso_score'] == 0 else row['fc'], axis=1)


#%%
# drug_counts = df['Drug'].value_counts()
# drug_counts_sorted = drug_counts.sort_values(ascending=False)
# drug_list = drug_counts_sorted.index.tolist()[:21]
# drug_df = df[df['Drug'].isin(drug_list)]
drug_df = df

#%%
# prot_counts = df['Gene'].value_counts()
# prot_counts_sorted = prot_counts.sort_values(ascending=False)
# prot_list = prot_counts_sorted.index.tolist()[:21]
# prot_df = df[df['Gene'].isin(prot_list)]
prot_df = df

#%%
with st.container():
    col1, col2, col3 = st.columns([1, 6, 1])
    with col1:
        pass
    
    with col2:
        # 中间列放置标题和作者
        st.markdown(
            "<h1 style='text-align: center; color: black; background-color: white;'>Our Title</h1>",
            unsafe_allow_html=True)

        st.markdown("""
        <p style='text-align: center; color: black; background-color: white;'>
        Our authors
        </p>
        """, unsafe_allow_html=True)

    with col3:
        # 右侧放置实验室logo
        # image = Image("C:/Users/17608/Desktop/MyAPP/figures/logo.jpg")
        st.image("images/logo.jpg", width=200)


#%%
st.subheader("Abstract")
st.write("""
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
""")

landscape, network1, network2 = st.tabs(["1️⃣ Scatter Plot", "2️⃣ Network Interaction (Drug)", "3️⃣ Network Interaction (Target)"])

#%%
with landscape:
    col1, col2, col3 = st.columns([1, 0.1, 1])
    selected_drug = col1.selectbox(":four_leaf_clover: Select a Drug", drug_df['Drug'].unique())
    selected_hscore = col3.slider(':herb: Set Hscore threshold', min_value=0.8, max_value=1.0, value=0.8, step=0.01)
    filtered_df = drug_df[(drug_df['Drug'] == selected_drug) & (drug_df['H-Score'] > selected_hscore)]

    # 绘制火山图
    col1, col2, col3 = st.columns([2, 5, 2])
    if not filtered_df.empty:
        fig = px.scatter(filtered_df, 
                            x='fc', 
                            y='log_pvalue', 
                            color='H-Score', 
                            labels={'x': 'fc', 'y': 'log_pvalue', 'target':'Gene'},
                            # color_continuous_scale=px.colors.diverging.RdBu,
                            color_continuous_scale="reds",
                            hover_data={'Gene': True}
                            )
        fig.update_traces(marker_size=12)
        fig.update_layout(
        margin=dict(l=0, r=0, t=2, b=0),  # 去掉左右和上下的边距
        autosize=True,
        width=1000,   # 设置图表的宽度
        height=450,
        # title_x=0.5,  # 设置标题居中
        xaxis_title_font=dict(size=20, weight='bold'),  # x轴标题字体大小
        yaxis_title_font=dict(size=20, weight='bold'),  # y轴标题字体大小
        xaxis_tickfont=dict(size=18),  # x轴刻度字体大小
        yaxis_tickfont=dict(size=18),)
        
        col2.plotly_chart(fig, use_container_width=True)
    else:
        col2.write('No results found.')


#%%
with network1:
    col11, col22, col33 = st.columns([1, 1, 1])
    selected_drug = col11.selectbox(":maple_leaf: Select a Drug", drug_df['Drug'].unique())
    layout= col22.selectbox(':fallen_leaf: Choose a network layout',('Random Layout','Spring Layout','Shell Layout','Kamada Kawai Layout'))
    # selected_hscore = col33.text_input('Set Hscore threshold')
    selected_hscore = col33.slider(':leaves: Set Hscore threshold', min_value=0.8, max_value=1.0, value=0.8, step=0.01)
    filtered_df = drug_df[(drug_df['Drug'] == selected_drug) & (drug_df['H-Score'] > selected_hscore)]

    # 绘制火山图
    col1, col2, col3 = st.columns([2, 5, 2])
    if not filtered_df.empty:
        
        # Create graph
        G = nx.Graph()

        # 添加节点和边 'Drug', 'Target', 'Hscore', 'logFC', 'logPvalue'
        
        for index, row in filtered_df.iterrows():
            if row['Drug'] not in G.nodes():
                G.add_node(row['Drug'], type='drug')  # Drug节点添加颜色和大小
            if row['Gene'] not in G.nodes():
                G.add_node(row['Gene'], type='target', score=row['H-Score'])  # Target节点为圆形
            if not G.has_edge(row['Drug'], row['Gene']):
                G.add_weighted_edges_from([(row['Drug'], row['Gene'], row['H-Score'])])
                
                # G.add_edge(row['Drug'], row['Targets'], weight=row['H-Score'])  
                # G.add_edge(row['Drug'], row['Targets'], weight=row['H-Score'])
        
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
                                                        title="H-Score",  # 颜色条的标题
                                                        titleside="top",  # 标题在右边
                                                        thickness=15,  # 设置颜色条的厚度
                                                        len=0.8,
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
                node_trace['marker']['symbol'] += tuple(['star'])   # Star shape for Drug
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
        margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
        autosize=False,
        width=500,  # 你可以设置图表的宽度
        height=500)
        col2.plotly_chart(fig, use_container_width=True)  

    else:
        col2.write(f"No data available for {selected_drug} with Hscore > {selected_hscore}.")


#%%
with network2:
    col111, col222, col333 = st.columns([1, 1, 1])
    selected_target = col111.selectbox(":seedling: Select a Target", prot_df['Gene'].unique())
    layout= col222.selectbox(':palm_tree: Choose a network layout',('Random Layout','Spring Layout','Shell Layout','Kamada Kawai Layout'))
    # selected_hscore = col33.text_input('Set Hscore threshold')
    selected_hscore = col333.slider(':chestnut: Set Hscore threshold', min_value=0.8, max_value=1.0, value=0.8, step=0.01)
    filtered_df = prot_df[(prot_df['Gene'] == selected_target) & (prot_df['H-Score'] > selected_hscore)]

    col1, col2, col3 = st.columns([2, 5, 2])
    if not filtered_df.empty:
        
        # Create graph
        G = nx.Graph()

        # 添加节点和边 'Drug', 'Target', 'Hscore', 'logFC', 'logPvalue'
        
        for index, row in filtered_df.iterrows():
            if row['Drug'] not in G.nodes():
                G.add_node(row['Drug'], type='drug', score=row['H-Score'])  # Drug节点添加颜色和大小
            if row['Gene'] not in G.nodes():
                G.add_node(row['Gene'], type='target')  # Target节点为圆形
            if not G.has_edge(row['Gene'], row['Drug']):
                G.add_weighted_edges_from([(row['Gene'], row['Drug'], row['H-Score'])])
                
        
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
                                            colorscale='reds',  # RdBu
                                            colorbar=dict(
                                                        title="H-Score",  # 颜色条的标题
                                                        titleside="top",  # 标题在右边
                                                        thickness=15,  # 设置颜色条的厚度
                                                        len=0.8,
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
            if node_type == 'target':
                node_trace['marker']['symbol'] += tuple(['circle'])   # Star shape for Drug
                node_trace['marker']['color'] += tuple(['#156082'])  # Red color for drug nodes
                node_trace['marker']['size'] += tuple([25])  # Fixed size for drug nodes
            elif node_type == 'drug':
                node_trace['marker']['symbol'] += tuple(['star'])
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
                                            # margin=dict(b=20,l=5,r=5,t=40),
                                            # margin=dict(b=5,l=5,r=5,t=5),
                                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

        fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
        autosize=False,
        width=500,  # 你可以设置图表的宽度
        height=500
        )
        col2.plotly_chart(fig, use_container_width=True)  

    else:
        col2.write(f"No data available for {selected_target} with Hscore > {selected_hscore}.")


#%%
st.divider()
st.subheader("Contact Us")
st.write("""🎈 Please feel free to contact us with any issues, comments, or questions.""")

st.text_area("**Message**")
st.text_input("**Your Email**")
if st.button("**Send**"):
    st.success("Message sent! We will get back to you soon.")


#  # 网页底部版权信息
# st.markdown(
#     '<div style="display: flex; justify-content: center; align-items: center; '
#     'height: 100px; background-color: #F2F2F2; color: black;">'
#     '<p style="margin: 0;">Copyright 1995-2024 Kanehisa Laboratories</p>'
#     '</div>',
#     unsafe_allow_html=True
# )