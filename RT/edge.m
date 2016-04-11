% The data structures required for the implementation of the 
% RT0 (Raviart Thomas Elements of the lowest order)

% Get the element data and the coordinates of the node from .dat files


% function [nodes2element, nodes2edge, noedges, edge2element,
%               interioredge] = edge(element,coordinates)

clear all
close all

format long

element = load('element.dat');
coordinate = load('coordinate.dat');

% nodes2element matrix : size(ndofs by ndofs)
nodes2element = sparse(size(coordinate,1),size(coordinate,1));
for j=1:size(element,1)
    nodes2element(element(j,:),element(j,[2 3 1]))= ...
        nodes2element(element(j,:),element(j,[2 3 1]))+j*eye(3,3);
end

%nodes2edge matrix : size(ndofs by ndofs)
B=nodes2element+nodes2element';
[I,J]=find(triu(B));
nodes2edge=sparse(I, J, 1:size(I,1),size(coordinate,1), size(coordinate,1));
nodes2edge=nodes2edge+nodes2edge';

%noedges vector (Number of edges vector)
noedges = size(I,1);

%edge2element : size(noedges,4)
edge2element=zeros(noedges,4);
for m=1:size(element,1)
    for k=1:3
        initial_node=element(m,k);
        end_node=element(m,rem(k,3)+1);
        p=nodes2edge(element(m,k),element(m,rem(k,3)+1));
        if edge2element(p,1)==0
            edge2element(p,:)=[initial_node ...
                              end_node ...
                              nodes2element(initial_node,end_node) ...
                              nodes2element(end_node,initial_node)];
        end
    end
end
