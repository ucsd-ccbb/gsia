##
## We need the igraph and Matrix libraries
##
.onLoad <- function(lib, pkg) {
   require(igraph)
   require(Matrix)
}



#' Correlation matrix G
#'
#' @param g undirected graph
#' @return numeric matrix G
#' @seealso ""
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28),ncol=2,byrow=TRUE))
#' species<-as.character(1:28)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=species)
#' set.seed(135)
#' plot(g,vertex.size=11)
#' G<-calculate_G(g)
calculate_G<-function(g) {
   vnames<-names(V(g))
   A<-as_adjacency_matrix(g)
   L<-diag(degree(g)+1)-A
   CC<-Matrix::tcrossprod(L,L)
   Ch<-Matrix::Cholesky(CC) #prepare Cholesky decomposition
   G<-as.matrix(Matrix::solve(Ch,L)) #inverse of L
   dG<-diag(G)
   G<-t(G/sqrt(dG))/sqrt(dG) #normalized G
   colnames(G)<-vnames
   rownames(G)<-vnames
   return(G)
}

#' Pairwise mutual information matrix I_{ij}
#'
#' @param G correlation matrix G
#' @return numeric matrix I_{ij}
#' @seealso ""
#' @export
#' @examples
#' Iij<-calculate_Iij(G)
calculate_Iij<-function(G) {
   u<-matrix(0,nrow=dim(G)[1],ncol=dim(G)[2])
   u[upper.tri(G)]<--1/2*log(1-G[upper.tri(G)]^2)
   u<-u+t(u)
   return(u)
}

#' Interaction information I(kA) of gene set A
#'
#' @param kA numeric vector of gene indices in set A
#' @param G correlation matrix G
#' @return numeric I(kA)
#' @seealso ""
#' @export
#' @examples
#' #here we are assuming that matrix G and graph g exist
#' A<-c(1,17,8) #vector of gene indices in set A
#' names(V(g))[A] # these are the corresponding gene names
#' tI(kA,G)
tI<-function(kA,G) { #interaction information
   d<-length(kA)
   u<-0
   if (d>1) u<- -1/2*determinant(G[kA,kA])[[1]][1]
   return(u)
}

#' Mutual information I(kA,kB) of gene sets A and B
#'
#' @param kA numeric vector of gene indices in set A
#' @param kB numeric vector of gene indices in set B
#' @param G correlation matrix G
#' @return numeric I(kA,kB)
#' @seealso ""
#' @export
#' @examples
#' indx<-1:length(V(g))
#' names(indx)<-as.character(V(g))
#' A<-c("1","17","8") #vector of gene names in set A
#' kA<-indx[A] #gene indices of set A
#' B<-c("1","11") #vector of gene names in set B
#' kB<-indx[B] #gene indices of set B
#' I(kA,kB,G)
I<-function(kA,kB,G) {
   u<-tI(base::union(kA,kB),G)-tI(setdiff(kA,kB),G)-tI(setdiff(kB,kA),G)
   return(u)
}

#' Variation of information VI(kA,kB) of a gene set A and B
#'
#' @param kA numeric vector of gene indices in set A
#' @param kB numeric vector of gene indices in set B
#' @param G correlation matrix G
#' @return numeric VI(kA,kB)
#' @seealso ""
#' @export
#' @examples
#' indx<-1:length(V(g))
#' names(indx)<-as.character(V(g))
#' A<-c("1","17","8") #vector of gene names in set A
#' kA<-indx[A] #gene indices in set A
#' B<-c("1","11") #vector of gene names in set B
#' kB<-indx[B] #gene indices in set B
#' VI(kA,kB,G)
VI<-function(kA,kB,G) {
   u<-tI(setdiff(kA,kB),G)+tI(setdiff(kB,kA),G)
   return(u)
}

#' Information quality ratio IQR(kA,kB) of a gene set A and B
#'
#' @param kA numeric vector of gene indices in set A
#' @param kB numeric vector of gene indices in set B
#' @param G correlation matrix G
#' @return numeric IQR(kA,kB)
#' @seealso ""
#' @export
#' @examples
#' indx<-1:length(V(g))
#' names(indx)<-as.character(V(g))
#' A<-c("1","17","8") #vector of gene names in set A
#' kA<-indx[A] #gene indices in set A
#' B<-c("1","11") #vector of gene names in set B
#' kB<-indx[B] #gene indices in set B
#' IQR(kA,kB,G)
IQR<-function(kA,kB,G) {
   u<-1-VI(kA,kB,G)/tI(base::union(kA,kB),G)
   return(u)
}
