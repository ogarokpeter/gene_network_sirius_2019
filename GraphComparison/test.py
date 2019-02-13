import graph2vec.trainer as gt

graph2vecc = gt.Graph2Vec(vector_dimensions=2)
for filename in ['edge_data.txt', 'edge_data2.txt', 'edge_data3.txt', 'edge_data4.txt', 'edge_data5.txt']:
    graph2vecc.parse_graph(filename, extend_paths=2)

    graph2vecc.fit(batch_size=1000, max_epochs=1000)

    print(filename, ":")
    print(graph2vecc.model.Win.get_value())
    print(graph2vecc.model.Wout.get_value())
