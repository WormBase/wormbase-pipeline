-- MySQL dump 10.13
--
-- Host: mcs2a    Database: mh6
-- ------------------------------------------------------
-- Server version	5.0.18-standard-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Not dumping tablespaces as no INFORMATION_SCHEMA.FILES table on this server
--

--
-- Table structure for table `CHROMOSOME_I`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_I` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `CHROMOSOME_II`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_II` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `CHROMOSOME_III`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_III` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `CHROMOSOME_IV`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_IV` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `CHROMOSOME_MtDNA`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_MtDNA` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `CHROMOSOME_V`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_V` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `CHROMOSOME_X`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `CHROMOSOME_X` (
  `tag_id` smallint(5) unsigned NOT NULL,
  `type_id` smallint(5) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL,
  KEY `start` (`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `chromo`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `chromo` (
  `tag_id` smallint(5) unsigned NOT NULL default '0',
  `type_id` smallint(5) unsigned NOT NULL default '0',
  `start` mediumint(8) unsigned NOT NULL default '0',
  `stop` mediumint(8) unsigned NOT NULL default '0',
  `frame` enum('.','0','1','2') NOT NULL default '.',
  `orientation` enum('+','-','.') NOT NULL default '.',
  `fluff` text NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `gff_tag`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `gff_tag` (
  `id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(50) NOT NULL default '',
  PRIMARY KEY  (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `gff_types`
--

SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `gff_types` (
  `id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(50) NOT NULL default '',
  PRIMARY KEY  (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2008-09-22 14:37:32
