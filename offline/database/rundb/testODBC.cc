#include <iostream>
#include <sstream>
#include <vector>
#include <typeinfo>
#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>

int main() {
    try {


      std::cout << "Testing ODBC." << std::endl;

      odbc::Connection* conn = odbc::DriverManager::getConnection("DSN=spinDB");
      std::cout << "Connected to spinDB." << std::endl;
      
      odbc::Statement* stmt = conn->createStatement();
      std::string query = "SELECT mbdvtx FROM spin WHERE runnumber = 45876";
      
      odbc::ResultSet* rs = stmt->executeQuery(query);

      while (rs->next()) {
	// ======== Using getString() ============ //
	std::cout << "getString method:" << std::endl;
	std::string mbdvtx = rs->getString("mbdvtx");
	std::cout << mbdvtx << std::endl;
	// ======================================= //

	// ======== Using getBytes() ============ //
	std::cout << "getBytes method:" << std::endl;
	odbc::Bytes bytes = rs->getBytes("mbdvtx");
	const signed char* data = bytes.getData();
	std::string strdata(reinterpret_cast<const char*>(data));
	std::cout << strdata << std::endl;
	// ======================================= //
      }

      delete rs;
      delete stmt;
      delete conn;
    } catch (odbc::SQLException& e) {
      std::cerr << "SQL Error: " << e.getMessage() << std::endl;
    }
    
    return 0;
}
