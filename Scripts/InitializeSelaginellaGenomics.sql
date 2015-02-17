
-- ---
-- Globals
-- ---

-- SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
-- SET FOREIGN_KEY_CHECKS=0;

-- ---
-- Table 'Sequences'
-- BLUELEAF sequences are from Trinity assemblies and include non-coding portions and chimeric sequences with multiple CDSs. All other sequences are CDS sequences only
-- ---

SET FOREIGN_KEY_CHECKS=0;

DROP TABLE IF EXISTS `Sequences`;
		
CREATE TABLE `Sequences` (
  `seqID` CHAR(25) NOT NULL,
  `sequence` MEDIUMTEXT NULL DEFAULT NULL,
  `species` CHAR(5) NULL DEFAULT NULL,
  `length` MEDIUMINT NULL DEFAULT NULL,
  PRIMARY KEY (`seqID`),
  KEY `species` (`species`),
  CONSTRAINT `sequences_ibfk_1` FOREIGN KEY (`species`) REFERENCES `Species` (`speciesID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- ---
-- Table 'Species'
-- 
-- ---

DROP TABLE IF EXISTS `Species`;
		
CREATE TABLE `Species` (
  `speciesID` CHAR(5) NOT NULL,
  `species` VARCHAR(255) NULL DEFAULT NULL,
  `genus` VARCHAR(255) NULL DEFAULT NULL,
  `source` VARCHAR(255) NULL DEFAULT NULL,
  PRIMARY KEY (`speciesID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- ---
-- Table 'CodingSequences'
-- Coding portions of transcripts, inferred from TransDecoder for non-1KP sequences. All other sequences are full length
-- ---

DROP TABLE IF EXISTS `CodingSequences`;
		
CREATE TABLE `CodingSequences` (
  `geneID` char(25) NOT NULL DEFAULT '' COMMENT 'identical to seqID for non-BLUELEAF sequences',
  `seqID` char(25) DEFAULT NULL,
  `species` char(5) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `strand` binary(1) DEFAULT NULL COMMENT '0 = -; 1= +',
  `start_codon` binary(1) DEFAULT NULL COMMENT '1 = 5 prime complete',
  `stop_codon` binary(1) DEFAULT NULL COMMENT '1 = 3 prime complete',
  PRIMARY KEY (`geneID`),
  KEY `seqID` (`seqID`),
  CONSTRAINT `codingsequences_ibfk_5` FOREIGN KEY (`seqID`) REFERENCES `Sequences` (`seqID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- ---
-- Table 'AthHomologs'
-- Based on Gramene biomart data
-- ---

DROP TABLE IF EXISTS `AthHomologs`;
		
CREATE TABLE `AthHomologs` (
  `id` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `seqID` CHAR(25) NULL DEFAULT NULL,
  `athID` CHAR(25) NULL DEFAULT NULL,
  `percen_identity` DECIMAL NULL DEFAULT NULL,
  `homology_type` CHAR(25) NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `seqID` (`seqID`),
  KEY `athID` (`athID`),
  CONSTRAINT `athhomologs_ibfk_2` FOREIGN KEY (`athID`) REFERENCES `Tair` (`athID`),
  CONSTRAINT `athhomologs_ibfk_1` FOREIGN KEY (`seqID`) REFERENCES `OrthoGroups` (`geneID`)

) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- ---
-- Table 'Counts'
-- Counts are inferred by running HTSeq on Tophat based mappings
-- ---

DROP TABLE IF EXISTS `Counts`;
		
CREATE TABLE `Counts` (
  `geneID` CHAR(25) NOT NULL DEFAULT '',
  `leaf1` INTEGER NULL DEFAULT NULL,
  `leaf2` INTEGER NULL DEFAULT NULL,
  `leaf3` INTEGER NULL DEFAULT NULL,
  `leaf4` INTEGER NULL DEFAULT NULL,
  PRIMARY KEY (`geneID`),
  CONSTRAINT `counts_ibfk_1` FOREIGN KEY (`geneID`) REFERENCES `CodingSequences` (`geneID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- ---
-- Table 'normalised'
-- Normalised expression data from DESeq
-- ---

DROP TABLE IF EXISTS `normalised`;
		
CREATE TABLE `normalised` (
  `geneID` CHAR(25) NOT NULL DEFAULT '',
  `leaf1` DECIMAL NULL DEFAULT NULL,
  `leaf2` DECIMAL NULL DEFAULT NULL,
  `leaf3` DECIMAL NULL DEFAULT NULL,
  `leaf4` DECIMAL NULL DEFAULT NULL,
  PRIMARY KEY (`geneID`),
  CONSTRAINT `normalised_ibfk_1` FOREIGN KEY (`geneID`) REFERENCES `CodingSequences` (`geneID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- ---
-- Table 'vsd'
-- variance stabalising transforamtion from DESeq
-- ---

DROP TABLE IF EXISTS `vsd`;
		
CREATE TABLE `vsd` (
  `geneID` CHAR(25) NOT NULL DEFAULT '',
  `leaf1` DECIMAL NULL DEFAULT NULL,
  `leaf2` DECIMAL NULL DEFAULT NULL,
  `leaf3` DECIMAL NULL DEFAULT NULL,
  `leaf4` DECIMAL NULL DEFAULT NULL,
  PRIMARY KEY (`geneID`),
  CONSTRAINT `vsd_ibfk_1` FOREIGN KEY (`geneID`) REFERENCES `CodingSequences` (`geneID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- ---
-- Table 'AthHomologs'
-- Based on Gramene biomart data
-- ---

DROP TABLE IF EXISTS `AthHomologs`;
		
CREATE TABLE `AthHomologs` (
  `id` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `seqID` CHAR(25) NULL DEFAULT NULL,
  `athID` CHAR(25) NULL DEFAULT NULL,
  `percen_identity` DECIMAL NULL DEFAULT NULL,
  `homology_type` CHAR(25) NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `seqID` (`seqID`),
  KEY `athID` (`athID`),
  CONSTRAINT `athhomologs_ibfk_2` FOREIGN KEY (`athID`) REFERENCES `Tair` (`athID`),
  CONSTRAINT `athhomologs_ibfk_1` FOREIGN KEY (`seqID`) REFERENCES `OrthoGroups` (`geneID`)

) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- ---
-- Table 'Tair'
-- 
-- ---

DROP TABLE IF EXISTS `Tair`;
		
CREATE TABLE `Tair` (
  `athID` CHAR(25) NOT NULL DEFAULT '',
  `description` MEDIUMTEXT NULL DEFAULT NULL,
  `symbol` VARCHAR(255) NULL DEFAULT NULL,
  PRIMARY KEY (`athID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



-- ---
-- Table 'ClusterInfo'
-- manually curated info about some of the OrthoGroups
-- ---

DROP TABLE IF EXISTS `ClusterInfo`;
		
CREATE TABLE `ClusterInfo` (
  `orthoID` INTEGER NULL AUTO_INCREMENT DEFAULT NULL,
  `function` MEDIUMTEXT NULL DEFAULT NULL,
  `notes` MEDIUMTEXT NULL DEFAULT NULL,
  `include` BINARY NULL DEFAULT NULL,
  `single-copy` BINARY NULL DEFAULT NULL,
  PRIMARY KEY (`orthoID`),
  CONSTRAINT `clusterinfo_ibfk_1` FOREIGN KEY (`orthoID`) REFERENCES `OrthoGroups` (`orthoID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- ---
-- Table 'OrthoGroups'
-- Ortholog groups are defined by orthoMCL and non-redundant sequences are identifed by PhyloTreePruner
-- ---

DROP TABLE IF EXISTS `OrthoGroups`;
		

CREATE TABLE `OrthoGroups` (
  `geneID` char(25) NOT NULL DEFAULT '',
  `orthoID` int(11) NOT NULL,
  `non_redundant` binary(1) DEFAULT NULL  COMMENT '1 = NR; 0 = redundant (based on PhyloTreePruner)',
  PRIMARY KEY (`geneID`),
  KEY `orthoID` (`orthoID`),
  CONSTRAINT `orthogroups_ibfk_1` FOREIGN KEY (`geneID`) REFERENCES `CodingSequences` (`geneID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



-- ---
-- Test Data
-- ---

-- INSERT INTO `Sequences` (`seqID`,`sequence`,`species`,`length`) VALUES
-- ('','','','');
-- INSERT INTO `Species` (`speciesID`,`species`,`genus`,`source`) VALUES
-- ('','','','');
-- INSERT INTO `OrthoGroups` (`geneID`,`orthoID`,`non_redundant`) VALUES
-- ('','','');
-- INSERT INTO `CodingSequences` (`geneID`,`seqID`,`species`,`start`,`end`,`strand`) VALUES
-- ('','','','','','');
-- INSERT INTO `AthHomologs` (`id`,`seqID`,`athID`,`percen_identity`,`homology_type`) VALUES
-- ('','','','','');
-- INSERT INTO `Tair` (`athID`,`description`,`symbol`) VALUES
-- ('','','');
-- INSERT INTO `Counts` (`geneID`,`leaf1`,`leaf2`,`leaf3`,`leaf4`) VALUES
-- ('','','','','');
-- INSERT INTO `normalised` (`geneID`,`leaf1`,`leaf2`,`leaf3`,`leaf4`) VALUES
-- ('','','','','');
-- INSERT INTO `vsd` (`geneID`,`leaf1`,`leaf2`,`leaf3`,`leaf4`) VALUES
-- ('','','','','');
-- INSERT INTO `ClusterInfo` (`orthoID`,`function`,`notes`,`include`,`single-copy`) VALUES
-- ('','','','','');

LOCK TABLES `Species` WRITE;
/*!40000 ALTER TABLE `Species` DISABLE KEYS */;

INSERT INTO `Species` (`speciesID`, `species`, `genus`, `source`)
VALUES
	('ABIJ','lepidophylla','Selaginella','1KP'),
	('ATH','thalliana','Arabidopsis','TAIR'),
	('EFJ','moellendorffii','Selaginella','ENSEMBL'),
	('JKAA','wallacei','Selaginella','1KP'),
	('KJYC','willdenowii','Selaginella','1KP'),
	('KRAUS','kraussiana','Selaginella','BLUELEAF'),
	('KUXM','selaginoides','Selaginella','1KP'),
	('LGDQ','apoda','Selaginella','1KP'),
	('MOEL','moellendorffii','Selaginella','BLUELEAF'),
	('SELMO','moellendorffii','Selaginella','JGI'),
	('UNC','uncinata','Selaginella','BLUELEAF'),
	('WILD','willdenowii','Selaginella','BLUELEAF'),
	('ZFGK','kraussiana','Selaginella','1KP'),
	('ZYCD','acanthonota','Selaginella','1KP'),
	('ZZOL','stauntoniana','Selaginella','1KP');

/*!40000 ALTER TABLE `Species` ENABLE KEYS */;
UNLOCK TABLES;

