#pragma once

#include <fstream>

inline void importBinary(std::ifstream& ifs, bool& val) { ifs.read((char*)(&val), sizeof(bool));  }
inline void exportBinary(std::ofstream& ofs, bool  val) { ofs.write((char*)(&val), sizeof(bool)); }

inline void importBinary(std::ifstream& ifs, char& val) { ifs.read((char*)(&val), sizeof(char));  }
inline void exportBinary(std::ofstream& ofs, char  val) { ofs.write((char*)(&val), sizeof(char)); }

inline void importBinary(std::ifstream& ifs, int& val) { ifs.read((char*)(&val), sizeof(int));  }
inline void exportBinary(std::ofstream& ofs, int  val) { ofs.write((char*)(&val), sizeof(int)); }

inline void importBinary(std::ifstream& ifs, unsigned int& val) { ifs.read((char*)(&val), sizeof(unsigned int));  }
inline void exportBinary(std::ofstream& ofs, unsigned int  val) { ofs.write((char*)(&val), sizeof(unsigned int)); }

//#ifdef _WIN64
//inline void importBinary(std::ifstream& ifs, std::size_t& val) { ifs.read((char*)(&val), sizeof(std::size_t));  }
//inline void exportBinary(std::ofstream& ofs, std::size_t  val) { ofs.write((char*)(&val), sizeof(std::size_t)); }
//#endif

inline void importBinary(std::ifstream& ifs, float& val) { ifs.read((char*)(&val), sizeof(float));  }
inline void exportBinary(std::ofstream& ofs, float  val) { ofs.write((char*)(&val), sizeof(float)); }

inline void importBinary(std::ifstream& ifs, double& val) { ifs.read((char*)(&val), sizeof(double));  }
inline void exportBinary(std::ofstream& ofs, double  val) { ofs.write((char*)(&val), sizeof(double)); }