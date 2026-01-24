#include "csvArray.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::csvArray::csvArray
(
    const word& csvName,
    const bool hasHeader,
    const scalar coef
)
:
    csvName_(csvName),
    hasHeader_(hasHeader),
    coef_(coef),
    array_(readCSV())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<scalarList> Foam::csvArray::readCSV()
{
    std::vector<std::vector<std::string>> data; // データを格納する2次元ベクトル

    std::ifstream file(csvName_);
    if (!file.is_open()) 
    {
        Info << "Failed to open the CSV file: " << csvName_ << endl;
        exit(FatalIOError);
    }

    std::string line;
    
    if (hasHeader_)
    {
    	std::getline(file, line); // ヘッダー行を読み飛ばす
    }
    
    while (std::getline(file, line))
    {
        std::vector<std::string> row; // 1行のデータを格納するベクトル
        std::istringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ',')) 
        {
            row.push_back(cell); // カンマで区切られたデータをベクトルに追加
        }

        data.push_back(row); // 行データを2次元ベクトルに追加
    }

    file.close(); // ファイルを閉じる
    
    List<scalarList> sList(data.size());
    label rowi =0;
    for (const auto& row : data)
    {
        scalarList rowList(row.size(),0);
        label celli = 0;
        for (const std::string& cell : row)
        {
            //std::cout << cell << "\t";
            rowList[celli] = coef_*stof(cell);
            celli++;
        }
        
        sList[rowi] = rowList;
        rowi++;
        //std::cout << std::endl;
    }
        
    return sList;
}



// ************************************************************************* //
