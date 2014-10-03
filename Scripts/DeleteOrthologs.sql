SET FOREIGN_KEY_CHECKS=0;
DROP TABLE IF EXISTS OrthoGroups;
CREATE TABLE `OrthoGroups` (
  `geneID` char(25) NOT NULL DEFAULT '',
  `orthoID` int(11) NOT NULL,
  `non_redundant` binary(1) DEFAULT NULL COMMENT '1 = NR; 0 = redundant (based on PhyloTreePruner)',
  PRIMARY KEY (`geneID`),
  KEY `orthoID` (`orthoID`),
  CONSTRAINT `orthogroups_ibfk_1` FOREIGN KEY (`geneID`) REFERENCES `CodingSequences` (`geneID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET FOREIGN_KEY_CHECKS=1;