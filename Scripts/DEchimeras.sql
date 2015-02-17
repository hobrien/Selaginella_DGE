SET FOREIGN_KEY_CHECKS=0;
DROP TABLE IF EXISTS `DEgenes`;
CREATE TABLE `DEgenes` (
  `key` int(11) NOT NULL AUTO_INCREMENT,
  `clusterID` char(25) NOT NULL DEFAULT '',
  `speciesID` char(5) DEFAULT NULL,
  `sample1` varchar(11) DEFAULT NULL,
  `sample2` varchar(11) DEFAULT NULL,
  `sample1Expr` float DEFAULT NULL,
  `sample2Expr` float DEFAULT NULL,
  `avgLogExpr` float DEFAULT NULL,
  `rLogFC` float DEFAULT NULL,
  `DGEclust_padj` float DEFAULT NULL,
  PRIMARY KEY (`key`),
  KEY `clusterID` (`clusterID`)
) ENGINE=InnoDB AUTO_INCREMENT=689042 DEFAULT CHARSET=latin1;
INSERT INTO DEgenes(clusterID, speciesID) 
    SELECT DISTINCT `clusterID`, `speciesID` FROM Expression 
        WHERE (rLogFC > 1 OR rLogFC < -1) AND DGEclust_padj < 0.01;
SET FOREIGN_KEY_CHECKS=1;
SELECT CodingSequences.seqID, CorsetGroups.clusterID, Count(CodingSequences.geneID) 
    FROM CodingSequences, CorsetGroups, DEgenes 
    WHERE CodingSequences.seqID = CorsetGroups.seqID 
        AND CorsetGroups.clusterID = DEgenes.clusterID 
        AND CodingSequences.species = DEgenes.speciesID 
    GROUP BY CodingSequences.seqID 
    HAVING Count(CodingSequences.geneID) > 1
    INTO OUTFILE '/tmp/chimeras.txt';