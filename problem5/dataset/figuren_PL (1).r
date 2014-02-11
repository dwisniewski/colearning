library(gplots)

to_dist <- function(K){
	n <- dim(K)[1]
		D <- K
	for(i in 1:n){
		for(j in 1:n){
			D[i,j] <- sqrt(K[i,i]+K[j,j]-2*K[i,j])
		}
	}
	D

}


KA <- read.table("~/Documents/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/y")
sum(KA>0)
#summary(t(KA))

KA <- array(data=as.matrix(KA),dim=c(127,38))

pdf("~/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/PL_figuren.pdf")

hist(KA, main="Histogram of the KA dataset")

ligands <- c('Roscovitine', 'Staurosporine', 'PKC-412', 'LY-333531', 'SB-202190', 'SU-14813', 'Sunitinib', 'CHIR-258/TKI-258', 'SB-203580', 'Flavopiridol', 'EKB-569', 'CP-724714', 'CP-690550', 'AMG-706', 'CHIR-265/RAF-265', 'MLN-8054', 'GW-2580', 'ABT-869', 'PI-103', 'GW-786034', 'JNJ-7706621', 'SB-431542', 'VX-745', 'MLN-518', 'BMS-387032/SNS-032', 'Sorafenib', 'Lapatinib', 'Erlotinib', 'BIRB-796', 'CI-1033', 'PTK-787', 'Gefitinib', 'AST-487', 'AZD-1152HQPA', 'Dasatinib', 'ZD-6474', 'Imatinib', 'VX-680/MK-0457')

proteins <- c('2G01', '1PMN', '3E7O', '2OU7', '2BUJ', '3COK', '3OCB', '2X39', '2W4O', '2CLQ', '2HW7', '2Z2W', '1U59', '3PP0', '3QQU', '3BBT', '2ITN', '2WEL', '2F57', '2V7O', '2VZ6', '3BHH', '2C47', '2IZS', '2CMW', '3LM0', '1OEC', '3E3B', '3RP0', '3GQI', '1XBC', '3H3C', '2O5K', '3EYG', '2CDZ', '2XIK', '2HY8', '3SOC', '2EB3', '2JBO', '2X4F', '2WQN', '2XNM', '2VWZ', '3DKO', '3A7H', '3H9F', '1SM2', '3G2F', '2Y7J', '2GDO', '3MY0', '3D7T', '1AD5', '1UWJ', '1U4D', '2ZV2', '2P2I', '3L8P', '2VAG', '3NR9', '2WU6', '1QPJ', '3OCS', '3GC7', '1CM8', '3FXX', '1MQB', '2X2L', '2X7F', '3QRK', '3FME', '3GVU', '3OMV', '3HNG', '3LCD', '3MTL', '3KN5', '2SRC', '2QLU', '3S95', '3LXN', '3GC8', '3CD3', '3BPR', '1WBP', '3A4O', '2X7G', '3LCS', '3C1X', '2CCH', '1UA2', '3RGF', '3BLQ', '2IJM', '3H9R', '2JC6', '2JAM', '3AQV', '2J7T', '3KRR', '3GT8', '3FAA', '2A19', '1GAG', '3AGL', '2DQ7', '1S9I', '3PP1', '3LL6', '3O50', '3ALN', '3BHY', '2CKE', '3EH9', '2IWI', '2XJ1', '2JFL', '1XJD', '3DTC', '2ZOQ', '2OJI', '2Z7R', '1T46', '1UWH', '3LXL', '3HRF')

PK <- read.table("~/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/Kernel_proteins", row.names=proteins)
hv <- heatmap.2(as.matrix(log(PK)), symm=TRUE,col= topo.colors(100), trace='none', main='Heatmap of the log. protein \n structure kernel', labRow=proteins, colCol=proteins)
PK <- to_dist(PK)

LK <-  read.table("~/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/Kernel_ligands", row.names=ligands)
hv <- heatmap.2(as.matrix((LK)), symm=TRUE,col= topo.colors(100), trace='none', main='Heatmap of the \n ligand kernel', labRow=ligands, labCol=ligands)
LK <- to_dist(LK)

plot(hclust(as.dist(LK)), main='Hierarchical clustering of \n the ligands')
plot(hclust(as.dist(PK)), main='Hierarchical clustering of \n the protein structures',cex=0.35)

Dock <- read.table('~/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/dockings')
Dock <- array(as.matrix(Dock),dim=c(127,38))

KA <- as.data.frame(KA, row.names=proteins)
names(KA) <- ligands
rowv <- as.dendrogram(hclust(dist(KA)))
colv <- as.dendrogram(hclust(dist(t(KA))))
plot(colv)
plot(rowv)

hv <- heatmap.2(as.matrix((KA)), col= topo.colors(100), trace='none',scale='row', main='Heatmap of the processed \n scaled KA dataset', Colv=colv, Rowv=rowv)

hv <- heatmap.2(as.matrix(log(log(KA+1)+1)), col= topo.colors(100), trace='none', main='Heatmap of the processed \n KA dataset', Colv=colv, Rowv=rowv)

Dock <- as.data.frame(Dock, row.names=proteins)
names(Dock) <- ligands

hv <- heatmap.2(as.matrix(Dock), col= topo.colors(100), trace='none', scale='row', main='Heatmap of the scaled \n docking data', Colv=colv, Rowv=rowv)


KA <-  read.table("~/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/y")
Dock <- read.table("~/Documents/Unief/Thesis/Protein-ligand_interaction/data_figuren/listD")
#plot(array(KA, dim=4826), array(Dock, dim =4826), main='KA data set vs. the docking results', xlab='KA',ylab='Docking')
dev.off()

