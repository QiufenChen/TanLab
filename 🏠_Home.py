import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx
# from PIL import Image


# Page configure
st.set_page_config(page_title="TanLab", layout="wide")
G = nx.random_geometric_graph(200, 0.125)

my_df = pd.read_csv("./Data/Dicoumarol.csv")
my_df['UniprotURL'] = ['https://www.uniprot.org/uniprotkb/{item}/entry' for item in my_df['Row.names']]
my_df['-log10P.Value'] = my_df['P.Value'].apply(lambda x: -np.log10(x))

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

    landscape, biomarker, network = st.tabs(["1️⃣ Landscape", "2️⃣ Biomarker", "🔥 Network"])
    with landscape:
        fig = px.scatter(
            my_df,
            x='logFC',
            y='-log10P.Value',
            color="Type",
            hover_data=['Protein.Ids', 'Genes', 'logFC', 'P.Value'],  # 在鼠标悬停时显示的内容
            labels={'UniprotID': 'Gene', 'log2FC': 'PValue'})

        fig.update_layout(
            margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
            autosize=False,
            width=600,  # 你可以设置图表的宽度
            height=600,
            title_x=0.5,  # 设置标题居中
            xaxis_title_font=dict(size=20),  # x轴标题字体大小
            yaxis_title_font=dict(size=20),  # y轴标题字体大小
            xaxis_tickfont=dict(size=18),  # x轴刻度字体大小
            yaxis_tickfont=dict(size=18)  # y轴刻度字体大小
        )

        # 添加两条竖线和一条横线
        fig.add_shape(
            type="line",
            x0=-0.1, x1=-0.1, y0=my_df['-log10P.Value'].min(), y1=my_df['-log10P.Value'].max(),
            line=dict(color="Red", width=2)
        )
        fig.add_shape(
            type="line",
            x0=0.1, x1=0.1, y0=my_df['-log10P.Value'].min(), y1=my_df['-log10P.Value'].max(),
            line=dict(color="Red", width=2)
        )
        fig.add_shape(
            type="line",
            x0=my_df['logFC'].min(), x1=my_df['logFC'].max(), y0=1, y1=1,
            line=dict(color="Blue", width=2)
        )
        # 在 Streamlit 中显示图表
        st.plotly_chart(fig, use_container_width=True)

    with biomarker:
        # reference ----> https://plotly.com/python/box-plots/
        swarm_df = px.data.tips()
        fig = px.strip(swarm_df, x='day', y='tip')
        fig.update_layout(
            margin=dict(l=0, r=0, t=0, b=0),  # 去掉左右和上下的边距
            autosize=False,
            width=500,  # 你可以设置图表的宽度
            height=400,
            title_x=0.5,  # 设置标题居中
            xaxis_title_font=dict(size=20),  # x轴标题字体大小
            yaxis_title_font=dict(size=20),  # y轴标题字体大小
            xaxis_tickfont=dict(size=18),  # x轴刻度字体大小
            yaxis_tickfont=dict(size=18)  # y轴刻度字体大小
        )
        st.plotly_chart(fig, use_container_width=True)

    with network:
        # Sample DataFrame Generation
        drugs = ['DrugA', 'DrugB', 'DrugC', 'DrugD', 'DrugE']
        data = []

        for drug in drugs:
            for i in range(1, 101):
                data.append([drug, f'Target{i}', random.uniform(0, 1)])

        df = pd.DataFrame(data, columns=['Drugs', 'Target', 'Hscore'])

        # Streamlit App
        st.title("Interactive Drug-Target Interaction Visualization")

        # Step 1: Drug Selection
        selected_drug = st.selectbox("Select a Drug", df['Drugs'].unique())

        # Step 2: Hscore threshold slider
        hscore_threshold = st.slider("Select Hscore threshold", min_value=0.0, max_value=1.0, value=0.8, step=0.01)

        # Step 3: Filter the DataFrame for the selected drug and Hscore threshold
        filtered_df = df[(df['Drugs'] == selected_drug) & (df['Hscore'] > hscore_threshold)]

        # Step 4: Display the filtered DataFrame
        st.write(f"Filtered Data for {selected_drug} with Hscore > {hscore_threshold}:")
        st.write(filtered_df)

        # Step 5: Create a graph (edges and nodes) using the filtered data
        if not filtered_df.empty:
            # Create graph
            G = nx.Graph()

            # Add nodes (drug and targets)
            G.add_node(selected_drug, type='drug', color='#C24A3A')  # Red pentagon for the drug
            for _, row in filtered_df.iterrows():
                target = row['Target']
                hscore = row['Hscore']
                G.add_node(target, type='target', color='#145975')  # Green circle for each target
                G.add_edge(selected_drug, target, weight=hscore)  # Create an edge with Hscore as weight

            # Step 6: Create a Plotly figure for network visualization
            pos = nx.spring_layout(G)

            # Create edge traces for each edge, with hover information for Hscore
            edge_x = []
            edge_y = []
            edge_text = []
            for edge in G.edges(data=True):
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
                edge_text.append(f'Hscore: {edge[2]["weight"]:.2f}')  # Hscore on hover

            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=2, color='blue'),
                hoverinfo='text',
                # text=edge_text,
                mode='lines'
            )

            # Create node traces for drugs and targets separately
            node_x = []
            node_y = []
            node_color = []
            node_text = []
            node_shape = []

            for node in G.nodes(data=True):
                x, y = pos[node[0]]
                node_x.append(x)
                node_y.append(y)
                node_color.append('red' if node[1]['type'] == 'drug' else 'green')
                node_text.append(node[0])  # Node name on hover
                node_shape.append('star' if node[1]['type'] == 'drug' else 'circle')

            node_trace = go.Scatter(
                x=node_x, y=node_y,
                mode='markers+text',
                hoverinfo='text',
                marker=dict(
                    showscale=False,
                    color=node_color,
                    size=20,
                    line_width=2,
                    symbol=['star' if shape == 'star' else 'circle' for shape in node_shape]  # Shape for drugs/targets
                ),
                # text=node_text
            )

            # Step 7: Layout and plot using Plotly
            fig = go.Figure(data=[edge_trace, node_trace],
                            layout=go.Layout(
                                title=f'Drug-Target Network for {selected_drug}',
                                titlefont_size=16,
                                showlegend=False,
                                hovermode='closest',
                                margin=dict(b=0, l=0, r=0, t=40),
                                xaxis=dict(showgrid=False, zeroline=False),
                                yaxis=dict(showgrid=False, zeroline=False))
                            )

            st.plotly_chart(fig)

        else:
            st.write(f"No data available for {selected_drug} with Hscore > {hscore_threshold}.")


    # --------------------------------------------------------------------------------------------------------------------------------------------------------
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
    st.subheader("Contact Us")
    st.write("""
    If you have any questions or need further information, feel free to reach out to us.
    """)

    st.text_area("Message")
    st.text_input("Your Email")
    if st.button("Send"):
        st.success("Message sent! We will get back to you soon.")
