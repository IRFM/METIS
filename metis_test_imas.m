% script for testing IMAS in METIS
diary off 
diary test_imasinmetis.log
diary on
zineb_path;
imasdb('test',getenv('USER'),'3');
metis4imas(1,1,'','auto_test');
exit