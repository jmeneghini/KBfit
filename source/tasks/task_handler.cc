#include "task_handler.h"
#include "stopwatch.h"
using namespace std;

// *************************************************************************

std::filesystem::path createKBOutputDirectory(const std::string& base_directory,
                                              const std::string& quantization_condition) {
  std::filesystem::path output_path;
  
  // Start with base directory if specified, this typically is the project name
  if (!base_directory.empty()) {
    output_path = std::filesystem::path(base_directory);
  }
  else {
    // current working directory
    output_path = std::filesystem::current_path();
  }

  // Add quantization condition subdirectory if specified
  if (!quantization_condition.empty()) {
    output_path = output_path / quantization_condition;
  }
  
  // Create the directory structure
  std::error_code ec;
  if (!std::filesystem::exists(output_path)) {
    if (!std::filesystem::create_directories(output_path, ec)) {
      if (ec) {
        throw(std::runtime_error("Error creating directory: " + ec.message()));
      }
    }
  }
  return output_path;
}

// set up the known tasks, create logfile, open stream for logging

TaskHandler::TaskHandler(XMLHandler& xmlin, int mpi_rank) : m_mpi_rank(mpi_rank) {
  if (xmlin.get_node_name() != "KBFit")
    throw(std::invalid_argument("Input file must have root tag <KBFit>"));
  string nowstr = get_date_time();

  if (xmlin.count_among_children("Initialize") != 1)
    throw(std::invalid_argument("There must be one child <Initialize> tag"));
  XMLHandler xmli(xmlin, "Initialize");

  if (xmli.count_among_children("MCSamplingInfo") != 1)
    throw(std::invalid_argument("There must be one <MCSamplingInfo> tag"));
  try {
    XMLHandler xmls(xmli, "MCSamplingInfo");
    MCSamplingInfo sampinfo(xmls);
    m_obs = new KBObsHandler(sampinfo);
  } catch (const std::exception& errmsg) {
    throw(std::invalid_argument("Failure reading sampling information"));
  }

  // get the user specified output dir
  if (xmli.count_among_children("OutputDirectory") == 1) {
    xmlread(xmli, "OutputDirectory", m_output_directory, "TaskHandler");
  } else
    m_output_directory = string(".");

  // Only rank 0 creates and writes to log files
  if (m_mpi_rank == 0) {
    string logfile;
    if (xmli.count_among_children("LogFile") == 1)
      xmlread(xmli, "LogFile", logfile, "TaskHandler");
    else
      logfile = string("kbfit_log_") + nowstr + ".xml";

    // put log file in the output directory
    
    filesystem::path logfile_path;
    if (logfile.find('/') == string::npos && logfile.find('\\') == string::npos) {
      filesystem::path project_dir = createKBOutputDirectory(m_output_directory);
      logfile_path = project_dir / logfile;
    } else {
      logfile_path = filesystem::path(logfile);
    }

    clog.open(logfile_path.c_str());
    if (!clog.is_open()) {
      cout << "Could not open log file " << logfile << " for output" << endl;
      throw(std::invalid_argument("Could not write to log file"));
    }
    clog << "<LogKBFit>" << endl;
  }

  string projname;
  if (xmli.count_among_children("ProjectName") == 1) {
    xmlread(xmli, "ProjectName", projname, "TaskHandler");
  } else
    projname = string("KBFit Project ") + nowstr;
  m_project_name = projname;
  
  // Only rank 0 writes to log
  if (m_mpi_rank == 0) {
    clog << " <ProjectName>" << projname << "</ProjectName>" << endl;
    clog << " <StartDateTime>" << nowstr << "</StartDateTime>" << endl;

    if (xmli.count_among_children("EchoXML") >= 1) {
      string input(xmlin.output());
      int pos = input.find("<KBFit>");
      input.erase(0, pos + 9);
      pos = input.find("</KBFit>");
      input.erase(pos, string::npos);
      clog << " <InputXML>" << endl;
      clog << input << "</InputXML>" << endl;
    }
  }

  m_task_map["DoPrint"] = &TaskHandler::doPrint;
  m_task_map["DoFit"] = &TaskHandler::doFit;
  m_task_map["DoSingleChannel"] = &TaskHandler::doSingleChannel;
}

// delete data, finish up logging, close log file

TaskHandler::~TaskHandler() {
  finish_log();
  m_task_map.clear();
  clear_task_data();
}

void TaskHandler::finish_log() {
  if (m_mpi_rank == 0) {
    string nowstr = get_date_time();
    clog << " <FinishDateTime>" << nowstr << "</FinishDateTime>" << endl;
    clog << "</KBFit>" << endl << endl;
    clog.close();
  }
}

void TaskHandler::do_batch_tasks(XMLHandler& xmlin) {
  XMLHandler xmlt(xmlin, "TaskSequence");
  list<XMLHandler> taskxml = xmlt.find_among_children("Task");
  int count = 0;
  
  if (m_mpi_rank == 0) {
    clog << endl
         << "<BeginTasks>****************************************</BeginTasks>"
         << endl;
  }
  
  for (list<XMLHandler>::iterator it = taskxml.begin(); it != taskxml.end();
       it++, count++) {
    if (m_mpi_rank == 0) {
      clog << endl << "<Task>" << endl;
      clog << " <Count>" << count << "</Count>" << endl;
    }
    
    XMLHandler xmlout;
    StopWatch rolex;
    rolex.start();
    do_task(*it, xmlout, count);
    rolex.stop();
    
    if (m_mpi_rank == 0) {
      clog << xmlout.output() << endl;
      clog << "<RunTimeInSeconds>" << rolex.getTimeInSeconds()
           << "</RunTimeInSeconds>" << endl;
      clog << "</Task>" << endl;
    }
  }
}

void TaskHandler::do_task(XMLHandler& xml_task, XMLHandler& xml_out,
                          int count) {
  try {
    XMLHandler xmlt(xml_task);
    if (xmlt.get_node_name() != "Task") {
      throw(std::invalid_argument("Input to do_task is not a Task tag"));
    }
    xmlt.seek_first_child();
    if (xmlt.get_node_name() != "Action") {
      throw(
          std::invalid_argument("Need Action tag as first child of Task tag"));
    }

    xmlt.seek_root();
    xml_child_assert(xmlt, "Action", "do_task");
    if (!xmlt.is_simple_element())
      throw(std::invalid_argument("Action tag is not simple XML element"));
    string task_action = xmlt.get_text_content();

    map<string, task_ptr>::iterator taskit = m_task_map.find(task_action);
    if (taskit != m_task_map.end()) {
      (this->*(taskit->second))(xml_task, xml_out, count);
    } // do the task!!
    else {
      throw(std::invalid_argument(string("Unknown task name: ") + task_action));
    }
  } // unknown task?
  catch (const std::exception& errmsg) {
    if (xml_out.empty())
      xml_out.set_root("Error", string(errmsg.what()));
    else
      xml_out.put_child("Error", string(errmsg.what()));
  }
}

string TaskHandler::get_date_time() {
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  string dt = asctime(timeinfo);
  for (unsigned int k = 0; k < dt.length(); k++)
    if (dt[k] == ' ')
      dt[k] = '_';
  return tidyString(dt);
}

void TaskHandler::clear_task_data() {
  for (map<string, TaskHandlerData*>::iterator it = m_task_data_map.begin();
       it != m_task_data_map.end(); it++)
    delete it->second;
  m_task_data_map.clear();
}

void TaskHandler::erase_task_data(const string& tdname) {
  string taskname = tidyName(tdname);
  map<string, TaskHandlerData*>::iterator it = m_task_data_map.find(taskname);
  if (it != m_task_data_map.end()) {
    delete it->second;
    m_task_data_map.erase(it);
  }
}

//  TaskHandlerData should be created with "new" so is persistent

void TaskHandler::insert_task_data(const std::string& tdname,
                                   TaskHandlerData* tdata) {
  string taskname = tidyName(tdname);
  if (taskname.empty()) {
    throw(
        std::invalid_argument("Cannot insert task data since name is invalid"));
  }
  map<string, TaskHandlerData*>::iterator it = m_task_data_map.find(taskname);
  if (it != m_task_data_map.end()) {
    throw(std::invalid_argument(
        "Cannot insert task data since name already in map"));
  }
  m_task_data_map.insert(make_pair(taskname, tdata));
}

//  returns a base pointer, which can be dynamic_cast(...)
//  Return null value if error.

TaskHandlerData* TaskHandler::get_task_data(const string& tdname) {
  string taskname = tidyName(tdname);
  if (taskname.empty())
    return 0;
  map<string, TaskHandlerData*>::iterator it = m_task_data_map.find(taskname);
  if (it == m_task_data_map.end())
    return 0;
  return it->second;
}

// *****************************************************************
/*
   //   <Task>
   //     <Action>ClearMemory</Action>
   //   </Task>

void TaskHandler::clearMemory(XMLHandler& xmltask, XMLHandler& xmlout, int
taskcount)
{
 m_obs->clearData();   // also clears all Samplings
 xmlout.set_root("ClearMemory","done");
}

   //   <Task>
   //     <Action>ClearSamplings</Action>
   //   </Task>

void TaskHandler::clearSamplings(XMLHandler& xmltask, XMLHandler& xmlout, int
taskcount)
{
 m_obs->clearSamplings();
 xmlout.set_root("ClearSamplings","done");
}

   //   <Task>
   //     <Action>EraseData</Action>
   //      <MCObservable>...</MCObservable>
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::eraseData(XMLHandler& xmltask, XMLHandler& xmlout, int
taskcount)
{
 xmlout.set_root("EraseData");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       MCObsInfo obskey(*tt);
       m_obs->eraseData(obskey);   // also clears all Samplings
       XMLHandler xmlb; obskey.output(xmlb);
       xmlout.put_child(xmlb);}
}

   //   <Task>
   //     <Action>EraseSamplings</Action>
   //      <MCObservable>...</MCObservable>
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::eraseSamplings(XMLHandler& xmltask, XMLHandler& xmlout, int
taskcount)
{
 xmlout.set_root("EraseSamplings");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       MCObsInfo obskey(*tt);
       m_obs->eraseSamplings(obskey);   // also clears all Samplings
       XMLHandler xmlb; obskey.output(xmlb);
       xmlout.put_child(xmlb);}
}


   //   <Task>
   //     <Action>ReadSamplingsFromFile</Action>
   //      <FileName>name_of_file</FileName>
   //      <MCObservable>...</MCObservable>   (these are optional)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::readSamplingsFromFile(XMLHandler &xmltask, XMLHandler& xmlout,
int taskcount)
{
 xmlout.set_root("ReadSamplingsFromFile");
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"TaskHandler");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 if (obskeys.empty())
    m_obs->readSamplingValuesFromFile(filename,xmlf);
 else
    m_obs->readSamplingValuesFromFile(obskeys,filename,xmlf);
 xmlout.put_child(xmlf);
}

   //   <Task>
   //     <Action>WriteSamplingsToFile</Action>   (uses samplings mode in
<MCSamplingInfo> tag)
   //      <FileName>name_of_file</FileName>
   //      <FileMode>overwrite</FileMode>   (optional)
   //      <MCObservable>...</MCObservable>   (these are needed)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::writeSamplingsToFile(XMLHandler &xmltask, XMLHandler& xmlout,
int taskcount)
{
 xmlout.set_root("WriteSamplingsToFile");
 string filename;
 xmlreadchild(xmltask,"FileName",filename,"TaskHandler");
 bool overwrite = false;  // protect mode
 if (xml_tag_count(xmltask,"FileMode")==1){
    string fmode;
    xmlread(xmltask,"FileMode",fmode,"FileListInfo");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") overwrite=true;}
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 m_obs->writeSamplingValuesToFile(obskeys,filename,xmlf,overwrite);
 xmlout.put_child(xmlf);
}



   //   <Task>
   //     <Action>ReadBinsFromFile</Action>
   //      <BinFileName>name_of_file</BinFileName>
   //      <MCObservable>...</MCObservable>   (these are optional)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::readBinsFromFile(XMLHandler &xmltask, XMLHandler& xmlout, int
taskcount)
{
 xmlout.set_root("ReadBinsFromFile");
 string filename;
 xmlreadchild(xmltask,"BinFileName",filename,"TaskHandler");
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 if (obskeys.empty())
    m_obs->readBinsFromFile(filename,xmlf);
 else
    m_obs->readBinsFromFile(obskeys,filename,xmlf);
 xmlout.put_child(xmlf);
}

   //   <Task>
   //     <Action>WriteBinsToFile</Action>
   //      <BinFileName>name_of_file</BinFileName>
   //      <FileMode>overwrite</FileMode>   (optional)
   //      <MCObservable>...</MCObservable>   (these are needed)
   //      <MCObservable>...</MCObservable>
   //   </Task>

void TaskHandler::writeBinsToFile(XMLHandler &xmltask, XMLHandler& xmlout, int
taskcount)
{
 xmlout.set_root("WriteBinsToFile");
 string filename;
 xmlreadchild(xmltask,"BinFileName",filename,"TaskHandler");
 bool overwrite = false;  // protect mode
 if (xml_tag_count(xmltask,"FileMode")==1){
    string fmode;
    xmlread(xmltask,"FileMode",fmode,"FileListInfo");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") overwrite=true;}
 list<XMLHandler> xmlh=xmltask.find("MCObservable");
 set<MCObsInfo> obskeys;
 for (list<XMLHandler>::iterator tt=xmlh.begin();tt!=xmlh.end();tt++){
       obskeys.insert(MCObsInfo(*tt));}
 XMLHandler xmlf;
 m_obs->writeBinsToFile(obskeys,filename,xmlf,overwrite);
 xmlout.put_child(xmlf);
}
*/
// ***************************************************************************************
