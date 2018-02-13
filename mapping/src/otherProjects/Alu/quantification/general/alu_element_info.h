// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALU_ELEMENT_INFO_H
#define ALU_ELEMENT_INFO_H

using namespace std;

class Alu_Element_Info
{
private:
	string id;
	string aluFamilyName;
	int chrNameInt;
	int startPos;
	int endPos;
	int length;

	// quantified results
	double readCount;
	double abundance;
public:
	Alu_Element_Info()
	{}

	void initiate(string& tmpId, string& tmpAluFamilyName, int tmpChrNameInt, int tmpStartPos, int tmpEndPos)
	{
		id = tmpId;
		aluFamilyName = tmpAluFamilyName;
		chrNameInt = tmpChrNameInt;
		startPos = tmpStartPos;
		endPos = tmpEndPos;
		length = endPos - startPos + 1;
		readCount = 0.0;
		abundance = 0.0;
	}

	void addReadCount(double toAddReadCount)
	{
		readCount = readCount + toAddReadCount;
	}

	string return_id()
	{
		return id;
	}

	int return_chrNameInt()
	{
		return chrNameInt;
	}

	int return_startPos()
	{
		return startPos;
	}

	int return_endPos()
	{
		return endPos;
	}

	int return_length()
	{
		return length;
	}

	double return_readCount()
	{
		return readCount;
	}

	double return_abundance()
	{
		return abundance;
	}
};
#endif