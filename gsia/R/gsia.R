##
## We need the igraph, mygene and Matrix libraries
##
.onLoad <- function(lib, pkg) {
   require(igraph)
   require(mygene)
   require(Matrix)
}


#' STRING PPI network annotated with Entrez gene IDs
#'
#' @param species taxonomy ID
#' @param score minimum interaction score
#' @return graph
#' @seealso ""
#' @export
#' @examples
#' #for human ppi graph with interaction score >= 700
#' g_STRING<-get_STRING() #for human ppi graph with interaction score >= 700
#'
#' #for mouse ppi graph with score >= 500
#' g_STRING_mouse<-get_STRING(species=10090,score=500) #for mouse ppi graph with score >= 500
get_STRING<-function(species=9606,score=700) {

   destfile<-paste(species,".protein.links.v10.5.txt.gz",sep="")
   download.file(paste("https://stringdb-static.org/download/protein.links.v10.5/",species,".protein.links.v10.5.txt.gz",sep=""),destfile=destfile)
   cat("reading",destfile,"...")
   X<-read.table(gzfile(destfile),sep=" ",header=TRUE)
   cat("done\n")

   good<-X$combined_score>=score
   Y<-X[good,1:2]
   p1<-sapply(strsplit(as.character(Y$protein1), "[.]"), "[", 2) #strip 9606.
   p2<-sapply(strsplit(as.character(Y$protein2), "[.]"), "[", 2)
   swap<-p1>p2 #to swap IDs
   p0<-p1
   p0[swap]<-p2[swap]
   p2[swap]<-p1[swap]
   df<-data.frame(p0,p2)
   df_unique<-unique(df)
   nrow(df_unique) #360341
   same<-as.character(df_unique[,1])==as.character(df_unique[,2])
   Z<-df_unique[!same,1:2]

   g<-graph_from_data_frame(Z,directed=FALSE)
   length(V(g)) #15154
   ENSP.names<-V(g)$name
   entrez.names<-queryMany(ENSP.names,scopes="ensembl.protein",fields="entrezgene",species=species)
   ENSP2entrez<-entrez.names[!duplicated(entrez.names$query),]
   good<-!is.na(entrez.names$entrezgene) & !duplicated(entrez.names$query)
   goodENSP2entrez<-entrez.names[good,]
   #sort genes by X_score and go with the better match
   o<-order(goodENSP2entrez$X_score,decreasing=TRUE)
   goodENSP2entrez_sorted<-goodENSP2entrez[o,]
   goodENSP2entrez_sorted_unique<-goodENSP2entrez_sorted[!duplicated(goodENSP2entrez_sorted$entrezgene),]
   #this is my official translation table with 14596 unique gene ids
   g0<-induced.subgraph(g,goodENSP2entrez_sorted_unique$query)
   g0<-simplify(g0) #remove loops
   mindx<-match(names(V(g0)),goodENSP2entrez_sorted_unique$query)
   V(g0)$name<-goodENSP2entrez_sorted_unique$entrezgene[mindx]
   #remove singletons
   single<-degree(g0)==0
   g0<-induced.subgraph(g0,!single)
   return(g0) #this is STRING translated into entrez gene
}

#' Calculate correlation matrix for a graph
#'
#' @param g undirected graph
#' @return numeric matrix G, Eq. (12) of Wallace et al. (2018)
#' @seealso ""
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
get_corr_matrix<-function(g) {
   vnames<-names(V(g))
   n<-length(V(g))
   t<-13*(n/14000)^3
   if (t>0.1) cat("this will take about ",format(t,digits=2,nsmall=1)," minutes...\n")
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
#' @return numeric matrix I_{ij}, Eq. (15) of Wallace et al. (2018).
#' @seealso ""
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' Iij<-get_Iij(G)
get_Iij<-function(G) {
   u<-G
   diag(u)<-0
   u<- -1/2*log(1-u^2)
   return(u)
}

#' Interaction information I(A) of gene set A
#'
#' @param A numeric or character vector of genes in set A
#' @param G correlation matrix G
#' @return numeric I(A), Eq. (13) of Wallace et al. (2018).
#' @seealso \code{\link{I}} for mutual information, \code{\link{VI}} for variation of information, and \code{\link{IQR}} for information quality ratio
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' A<-c("2","4","6","30")
#' tI(A,G)
tI<-function(A,G) { #interaction information
   d<-length(A)
   u<-0
   if (d>1) u<- -1/2*determinant(G[A,A])[[1]][1]
   return(u)
}

#' Mutual information I(A,B) of gene sets A and B
#'
#' @param A numeric or character vector of genes in set A
#' @param B numeric or character vector of genes in set B
#' @param G correlation matrix G
#' @return numeric I(A,B)
#' @seealso \code{\link{tI}} for interaction information, \code{\link{VI}} for variation of information, and \code{\link{IQR}} for information quality ratio
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' A<-c("2","4","6","30")
#' B<-c("21","23","29")
#' I(A,B,G)
I<-function(A,B,G) {
   u<-tI(base::union(A,B),G)-tI(setdiff(A,B),G)-tI(setdiff(B,A),G)
   return(u)
}

#' Variation of information VI(A,B) of gene sets A and B
#'
#' @param A numeric or character vector of genes in set A
#' @param B numeric or character vector of genes in set B
#' @param G correlation matrix G
#' @return numeric VI(A,B)
#' @seealso \code{\link{tI}} for interaction information, \code{\link{I}} for mutual information, and \code{\link{IQR}} for information quality ratio
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' A<-c("2","4","6","30")
#' B<-c("21","23","29")
#' VI(A,B,G)
VI<-function(A,B,G) {
   u<-tI(setdiff(A,B),G)+tI(setdiff(B,A),G)
   return(u)
}

#' Information quality ratio IQR(A,B) of gene sets A and B
#'
#' @param A numeric or character vector of genes in set A
#' @param B numeric or character vector of genes in set B
#' @param G correlation matrix G
#' @return numeric IQR(A,B)
#' @seealso \code{\link{tI}} for interaction information, \code{\link{I}} for mutual information, and \code{\link{VI}} for variation of information
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' A<-c("2","4","6","30")
#' B<-c("21","23","29")
#' IQR(A,B,G)
IQR<-function(A,B,G) {
   u<-1-VI(A,B,G)/tI(base::union(A,B),G)
   return(u)
}

#' Graph of expressed.genes as verices with interactions inherited from graph g
#' 
#' @param g graph
#' @param expressed.genes character vector of expressed gene IDs
#' @return graph containing all expressed.genes as vertices with interactions inherited from graph g. If any expressed.genes are not found in g, they will become unlinked vertices 
#' @seealso ""
#' @export
#' @examples
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' expressed.genes<-c("2","4","6","30","40")
#' g.expressed<-inherit_interactions(g,expressed.genes)
#' plot(g.expressed,vertex.size=11)
inherit_interactions<-function(g,expressed.genes) {
   vnames<-names(V(g))
   expressed.genes.in<-intersect(expressed.genes,vnames)
   expressed.genes.out<-setdiff(expressed.genes,vnames)
   g.expressed.in<-induced.subgraph(g,expressed.genes.in)
   g.expressed<-add.vertices(g.expressed.in,length(expressed.genes.out),name=expressed.genes.out) #works even if length==0
   return(g.expressed)
}

#' p-value of tI(A)
#'
#' @param A character vector of genes in set A
#' @param G correlation matrix G
#' @param null.domain gene set from which null gene sets are drawn (default colnames(G))
#' @param null.prob genes will be drawn with probability proportional to these weights (default is uniform probabilities).
#' @param nsamples number of null sets to be generated (default 1000)
#' @param plot logical. Should the null distribution be plotted (default TRUE)
#' @return The probability that tI(R) >= tI(A), where R is a null gene set of size |A|
#' drawn randomly from null.domain.
#' Each gene is drawn without replacement with a probability proportional to null.prob.
#' A total of nsamples will be drawn.
#' 
#' @seealso \code{\link{get_p_I}}
#' @export
#' @examples 
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' A<-c("2","4","6","30")
#' p<-get_p_tI(A,G)
get_p_tI<-function(A,G,null.domain=colnames(G),null.prob=rep(1,length(null.domain)),nsamples=1e3,plot=TRUE) {
   vnames<-colnames(G)
   A0<-intersect(A,vnames)
   Ai<-setdiff(A,vnames)
   if (length(Ai)>0) cat(length(Ai),"genes from A are not in the graph and will be ignored.\n")
   n<-length(A0)
   tI_null<-rep(0,nsamples)
   tIA<-tI(A0,G)
   indx_null<-match(null.domain,colnames(G))
   for (i in 1:nsamples) {
      A_null<-sample(indx_null,size=n,prob=null.prob) #indices, not names
      tI_null[i]<-tI(A_null,G)
   }
   p<-mean(tI_null>=tIA)
   if (plot) {
      r<-range(tI_null)
      hist(tI_null,breaks=100,xlim=c(min(tIA,r[1]),max(tIA,r[2])),xlab="null tI(A)",main="")
      rug(tI_null)
      rug(tIA,lwd=2,col="red")
      box()
   }
   return(p)
}

#' p-value of I(A,B)
#'
#' @param A character vector of genes in set A
#' @param B character vector of genes in set B
#' @param G correlation matrix G
#' @param null.domain gene set from which null gene sets are drawn (default colnames(G))
#' @param null.prob genes will be drawn with probability proportional to these weights (default is uniform probabilities)
#' @param nsamples number of null sets to be generated (default 1000)
#' @param plot logical. Should the null distribution be plotted (default TRUE)
#' @return The probability that I(A,R) >= I(A,B), where R is a null gene set of size |B|
#' drawn randomly from null.domain.
#' Each gene is drawn without replacement with a probability proportional to null.prob.
#' A total of nsamples will be drawn.
#' 
#' @seealso \code{\link{get_p_tI}}
#' @export
#' @examples 
#' df0<-as.data.frame(matrix(c(1,2,2,3,3,4,4,5,5,6,4,7,5,7,3,7,7,6,3,8,8,9,6,10,1,11,2,7,11,13,11,14,11,15,11,16,11,17,11,18,11,19,11,20,19,20,20,21,2,22,6,2,17,8,9,23,23,24,12,24,23,25,25,26,25,27,27,28,29,30),ncol=2,byrow=TRUE))
#' vnames<-as.character(1:30)
#' g<-graph_from_data_frame(df0,directed=FALSE,vertices=vnames)
#' set.seed(1357)
#' plot(g,vertex.size=11)
#' G<-get_corr_matrix(g)
#' A<-c("2","4","6","30")
#' B<-c("7","3","8","28","29")
#' p<-get_p_I(A,B,G)
get_p_I<-function(A,B,G,null.domain=colnames(G),null.prob=rep(1,length(null.domain)),nsamples=1e3,plot=TRUE) {
   vnames<-colnames(G)
   A0<-intersect(A,vnames)
   Ai<-setdiff(A,vnames)
   if (length(Ai)>0) cat(length(Ai),"genes from A are not in the graph and will be ignored.\n")
   B0<-intersect(B,vnames)
   Bi<-setdiff(B,vnames)
   if (length(Bi)>0) cat(length(Bi),"genes from B are not in the graph and will be ignored.\n")
   n<-length(B0)
   I_null<-rep(0,nsamples)
   IAB<-I(A0,B0,G)
   indx_null<-match(null.domain,colnames(G))
   indx_A<-match(A0,colnames(G))
   for (i in 1:nsamples) {
      B_null<-sample(indx_null,size=n,prob=null.prob) #indices, not names
      I_null[i]<-I(indx_A,B_null,G)
   }
   p<-mean(I_null>=IAB)
   if (plot) {
      r<-range(I_null)
      hist(I_null,breaks=100,xlim=c(min(IAB,r[1]),max(IAB,r[2])),xlab="null I(A,B)",main="")
      rug(I_null)
      rug(IAB,lwd=2,col="red")
      box()
   }
   return(p)
}
