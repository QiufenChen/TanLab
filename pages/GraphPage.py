#!/usr/bin/env python
# Author  : KerryChen
# File    : data.py
# Time    : 2024/9/11 17:20
import streamlit as st
import plotly.graph_objects as go
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# # Streamlit App
# st.title("Drug-Target Interaction Visualization")

# # 创建5个药物的列表
# drugs = ['DrugA', 'DrugB', 'DrugC', 'DrugD', 'DrugE']

# # 每个药物对应100个target
# targets = [f'Target{i}' for i in range(1, 101)]

# # 生成DataFrame
# data = []
# for drug in drugs:
#     for target in targets:
#         hscore = np.random.rand()
#         data.append([drug, target, hscore])

# # 创建DataFrame
# df = pd.DataFrame(data, columns=['Drugs', 'Target', 'HScore'])

# # Step 1: Drug Selection
# select_drug = st.selectbox("Select a Drug", df['Drugs'].unique())

# # Step 2: Hscore threshold slider
# hscore_threshold = st.slider("Select HScore threshold", min_value=0.0, max_value=1.0, value=0.8, step=0.1)


# # Step 2: Filter the DataFrame for the selected drug
# filtered_df = df[(df['Drugs'] == select_drug) & (df['HScore'] > hscore_threshold)]
# # Step 3: Display the filtered DataFrame
# st.write(f"Filtered Data for {select_drug} with HScore > {hscore_threshold}:")
# st.write(filtered_df)

# # Step 5: Create a graph (edges and nodes) using the filtered data
# if not filtered_df.empty:
#     # Create graph
#     G = nx.Graph()

#     # Add nodes (drug and targets)
#     G.add_node(select_drug, color='red', node_shape='p')  # Red pentagon for the drug
#     for _, row in filtered_df.iterrows():
#         target = row['Target']
#         hscore = row['HScore']
#         G.add_node(target)  # Green circle for each target
#         G.add_edge(select_drug, target, weight=hscore)     # Create an edge with Hscore as weight

#     # Step 6: Draw the graph with edge colors based on Hscore
#     st.subheader(f"Network for {select_drug} with HScore > {hscore_threshold}")

#     # Define edge colors based on Hscore
#     edges = G.edges(data=True)
#     edge_colors = [data['weight'] for _, _, data in edges]

#     # Create a figure and axis for plotting
#     fig, ax = plt.subplots()

#     # Define positions for the nodes
#     pos = nx.spring_layout(G)

#     # Draw the drug node (red pentagon)
#     drug_nodes = [select_drug]
#     nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color='#C24A3A', node_size=400, node_shape='*', ax=ax)

#     # Draw the target nodes (green circles)
#     target_nodes = [node for node in G.nodes if node != select_drug]
#     nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_color='#145975', node_size=300, node_shape='o', ax=ax)

#     # # Draw the labels
#     # nx.draw_networkx_labels(G, pos, font_size=10, ax=ax)

#     # Map edge colors based on the Hscore (using a colormap)
#     cmap = plt.cm.get_cmap('PuBu')
#     norm = plt.Normalize(vmin=min(edge_colors), vmax=max(edge_colors))

#     # Draw edges with varying colors based on Hscore
#     nx.draw_networkx_edges(G, pos, edge_color=edge_colors, edge_cmap=cmap, edge_vmin=min(edge_colors),
#                            edge_vmax=max(edge_colors), width=2, ax=ax)

#     # Add colorbar to indicate Hscore values
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm.set_array(edge_colors)
#     fig.colorbar(sm, ax=ax, label='HScore')

#     # Display the network graph in Streamlit
#     st.pyplot(fig)

# else:
#     st.write(f"No data available for {select_drug} with Hscore > {hscore_threshold}.")


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
