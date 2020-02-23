# coding:utf-8
# 2019.12.16
# 按照k为3分的part中，会出现有三个元素的part将其分开
# 按照一系列程序进行切分后，会得到一些小的part，根据他与其他part的overlap和相连情况拉倒同一个part中，
import sys
import networkx as nx
import time
from networkx.algorithms.community import k_clique_communities


def get_point_dict(input_line, line_dict1):
    for two_num in input_line:
        if int(two_num[0]) in line_dict1.keys():
            line_dict1[int(two_num[0])].append(int(two_num[1]))
        else:
            line_dict1[int(two_num[0])] = []
            line_dict1[int(two_num[0])].append(int(two_num[1]))
        if int(two_num[1]) in line_dict1.keys():
            line_dict1[int(two_num[1])].append(int(two_num[0]))
        else:
            line_dict1[int(two_num[1])] = []
            line_dict1[int(two_num[1])].append(int(two_num[0]))


def frozenset_list(set_list):
    index = 0
    set_list1 = []
    for set1 in set_list:
        set_list1.append(list(set1))
        index += 1
    return set_list1


def make_graph(line, k_value,set_num):
    # 创建一个无向、无权图
    # 定义图
    Graph = nx.Graph()
    # 获取边列表edges_list
    edges_list = []
    for edge in line:
        edge_list = edge
        edges_list.append((int(edge_list[0]), int(edge_list[1])))
    # 为图增加边
    Graph.add_edges_from(edges_list)
    part_list = list(k_clique_communities(Graph, k_value))
    part_list = frozenset_list(part_list)
    # print ('k为3切分：' + str(part_list))
    # 对part_list中part基因个数超标的行进行以k为四进行重新检验
    loop_index = 0
    loop_num = len(part_list)
    for part in part_list:
        loop_index += 1
        if loop_index > loop_num:
            break
        if len(part) > set_num:
            Graph2 = nx.Graph()
            edges_list = []
            for edge in line:
                if (int(edge[0]) in part) and (int(edge[1]) in part):
                    edge_list = edge
                    edges_list.append((int(edge_list[0]), int(edge_list[1])))
            Graph2.add_edges_from(edges_list)
            check_list = list(k_clique_communities(Graph2, 4))
            # print ('k为4切分' + str(check_list))
            # 对check后的列表进行一个判断，
            if len(check_list) == 1:
                continue
            else:
                part_list.remove(part)
                part_list += check_list
    # 将切分得到的part，进行重新的check一遍，将有overla和连线的part连在一起
    line_dict = {}
    get_point_dict(line, line_dict)
    check_list_sorted = sorted(part_list, key=lambda g: len(g), reverse=False)
    check_list_sorted = frozenset_list(check_list_sorted)
    check_index = -1
    del_check_list = []
    # 将part从少到多与其他part进行比较
    for check_part in check_list_sorted[:-1]:
        check_index += 1
        suit_num_list = []
        check_part1_index = check_index
        for check_part1 in check_list_sorted[check_index + 1:]:
            check_part1_index += 1
            suit_num = 0
            suit_list = []
            # 看下这个part里有多少点与被比part有关系，同时记录下被比part里有多少点被比上
            # 取两者的最少数
            for check_point in check_part:
                if check_point in check_part1:
                    suit_list.append(check_point)
                    suit_num += 1
                    continue
                for check_connect_point in line_dict[check_point]:
                    if check_connect_point in check_part1:
                        suit_list.append(check_connect_point)
                        suit_num += 1
                        break
            suit_num_list.append([min(suit_num, len(set(suit_list))), check_part1_index])
        suit_num_list_sorted = sorted(suit_num_list, key=lambda p: p[0], reverse=True)
        # 如果有关系的点占据这个part的一半的话，就认为他俩可以合并，
        if suit_num_list_sorted[0][0] * 2 >= len(check_part):
            check_list_sorted[suit_num_list_sorted[0][1]] += check_part
            check_list_sorted[suit_num_list_sorted[0][1]] = list(
                set(check_list_sorted[suit_num_list_sorted[0][1]]))
            del_check_list.append(check_index)
    del_num = -1
    for del_index in del_check_list:
        del_num += 1
        del check_list_sorted[del_index - del_num]
    # print ('合并后:' + str(check_list_sorted))
    # 分别遍历这个社区列表的所有part，再把part中可能因为达不到最小子图而踢掉的基因拉进来，并判断下他还有没出现在别的part里
    # 如果只是一个part或者就没有part的话，就把所有的基因都放进这个part
    part_list=check_list_sorted
    write_list = []
    if len(part_list) == 1:
        for two_num in line:
            for i in two_num:
                if i == '':
                    continue
                if int(i) in part_list[0]:
                    continue
                else:
                    part_list[0] = list(part_list[0])
                    part_list[0].append(int(i))
        write_list = part_list
    elif len(part_list) == 0:
        part_list.append([])
        for two_num in line:
            for i in two_num:
                if i in part_list[0]:
                    continue
                else:
                    part_list[0].append(i)
        write_list = part_list
    else:
        #  遍历这行将所有的基因序号和他的共线性基因序号写入字典
        line_dict = {}
        total_part = []
        for part in part_list:
            total_part += part
        get_point_dict(line, line_dict)
        # 拿到所有写进part的基因序号列表
        for part in part_list:
            part = list(part)
            part_add = []
            other_part = total_part.copy()
            for n in part:
                other_part.remove(n)
            for i in part:
                total_j_list = []
                if i in other_part:
                    continue
                for j in line_dict[i]:
                    if (j not in total_part) and (j not in total_j_list):
                        total_j_list.append(j)
                for m in total_j_list:
                    total_m_list = [m]
                    for x in total_m_list:
                        if x in part:
                            continue
                        elif x in total_part:
                            total_m_list = []
                            break
                        else:
                            for h in line_dict[x]:
                                if h not in total_m_list:
                                    if h not in part:
                                        total_m_list.append(h)
                    part_add += total_m_list
            part += part_add
            write_list.append(part)
    write_index = 0
    for write_part in write_list:
        write_list[write_index] = list(set(write_part))
        write_index += 1
    return write_list


if __name__ == '__main__':
    make_graph()