// Fill out your copyright notice in the Description page of Project Settings.


#include "GenerationFun.h"


float  UGenerationFun::erosionScore(const TArray<FVector> &coordVert, TArray<float> &matrixS)
{
	matrixS.Empty();
	int sizeMatrix = sqrt(coordVert.Num());
	int i, j;
	int cur = 0;
	float maxDif; // max differance
	float summ = 0;
	for (i = 0; i < sizeMatrix; i++)
	{
		for (j = 0; j < sizeMatrix; j++)
		{
			cur = i * sizeMatrix + j;// current cell
			//find max differance
			maxDif = 0;
			if ((i-1 >= 0) && (maxDif < abs(coordVert[cur].Z - coordVert[cur - sizeMatrix].Z)))// h_(i-1,j)
				maxDif = abs(coordVert[cur].Z - coordVert[cur - sizeMatrix].Z);
			if ((i+1 < sizeMatrix) && (maxDif < abs(coordVert[cur].Z - coordVert[cur + sizeMatrix].Z)))// h_(i+1,j)
				maxDif = abs(coordVert[cur].Z - coordVert[cur + sizeMatrix].Z);
			if ((j-1 >= 0) && (maxDif < abs(coordVert[cur].Z - coordVert[cur - 1].Z)))// h_(i,j-1)
				maxDif = abs(coordVert[cur].Z - coordVert[cur - 1].Z);
			if ((j+1 < sizeMatrix) && (maxDif < abs(coordVert[cur].Z - coordVert[cur + 1].Z)))// h_(i,j+1)
				maxDif = abs(coordVert[cur].Z - coordVert[cur + 1].Z);
			matrixS.Add(maxDif);
			summ += maxDif;
			
		}
	}
	float amountS = summ / pow(sizeMatrix,2);
	float dS = 0;
	for (i = 0; i < matrixS.Num(); i++)
	{
		dS += pow(matrixS[i] - amountS , 2);
	}
	dS = sqrt(dS / pow(sizeMatrix, 2));

	float eScore = dS / amountS;

	return eScore;
}



void UGenerationFun::SaveDataToFile(const TArray<FString>& Data, FString Name)
{
	FString GameDir = FPaths::ProjectContentDir();
	FString CompleteFilePath = GameDir + Name;
	FString ResultString;
	for (int i = 0; i < Data.Num(); i++)
		ResultString += Data[i];
	FFileHelper::SaveStringToFile(ResultString, *CompleteFilePath);
}

FString UGenerationFun::NextString(FString Data) {
	return(Data + "\n");
}

TArray<FVector> UGenerationFun::changeVertexArrayCoord(const TArray<FVector>& arrayVertex,float coefXY, float coefZ) {
	TArray<FVector> NewArrayVertex;
	int size = sqrt(arrayVertex.Num());
	FVector CurrentPoint;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			CurrentPoint.X = i * coefXY;
			CurrentPoint.Y = j * coefXY;
			CurrentPoint.Z = arrayVertex[i*size+j].Z * coefZ;
			NewArrayVertex.Add(CurrentPoint);
		}
	}
	return NewArrayVertex;
}

TArray<FVector> UGenerationFun::CreateVertexArrayComponent(const TArray<FVector>& arrayVertex, int indexComponent, int width) {
	TArray<FVector> ArrayVertexComponent;
	int i, j;
	int size = sqrt(arrayVertex.Num());
	int idInColl = size / width;
	int currentIdColl = indexComponent % idInColl;
	int currentIdRow = indexComponent / idInColl;

	for (i = 0; i < width; i++) {
		for (j = 0; j < width; j++) {
			ArrayVertexComponent.Add(arrayVertex[(currentIdRow*(width-1) + i)*size + currentIdColl*(width-1) + j]);
		}
	}
	return ArrayVertexComponent;
}

TArray<FVector> UGenerationFun::CreateTrees(const TArray<FVector>& coordVert) {
	
	TArray<float> matrixS;
	matrixS.Empty();
	int sizeMatrix = sqrt(coordVert.Num());
	int i, j;
	int cur = 0;
	float maxDif; // max differance
	float summ = 0;
	for (i = 0; i < sizeMatrix; i++)
	{
		for (j = 0; j < sizeMatrix; j++)
		{
			cur = i * sizeMatrix + j;// current cell
			//find max differance
			maxDif = 0;
			if ((i - 1 > 0) && (maxDif < abs(coordVert[cur].Z - coordVert[cur - sizeMatrix].Z)))// h_(i-1,j)
				maxDif = abs(coordVert[cur].Z - coordVert[cur - sizeMatrix].Z);
			if ((i + 1 < sizeMatrix) && (maxDif < abs(coordVert[cur].Z - coordVert[cur + sizeMatrix].Z)))// h_(i+1,j)
				maxDif = abs(coordVert[cur].Z - coordVert[cur + sizeMatrix].Z);
			if ((j - 1 > 0) && (maxDif < abs(coordVert[cur].Z - coordVert[cur - 1].Z)))// h_(i,j-1)
				maxDif = abs(coordVert[cur].Z - coordVert[cur - 1].Z);
			if ((j + 1 < sizeMatrix) && (maxDif < abs(coordVert[cur].Z - coordVert[cur + 1].Z)))// h_(i,j+1)
				maxDif = abs(coordVert[cur].Z - coordVert[cur + 1].Z);
			matrixS.Add(maxDif);
			summ += maxDif;
		}
	}

	float amountS = summ / pow(sizeMatrix, 2);
	TArray<FVector> spawnTree;
	for (i = 0; i < matrixS.Num(); i++) {
		if (matrixS[i] < 5) {
			spawnTree.Add(coordVert[i]);
		}
	}
	return spawnTree;
}

TArray<FVector> UGenerationFun::NormalizeCoordVector(TArray<FVector> coordVert) {
	int i;
	float max = 0;
	float coef;
	for (i = 0; i < coordVert.Num(); i++)
		if (coordVert[i].Z > max)
			max = coordVert[i].Z;
	//UE_LOG(LogTemp, Warning, TEXT("max = %f"), max);
	coef = 1 / max;
	//UE_LOG(LogTemp, Warning, TEXT("coef = %f"), coef);
	for (i = 0; i < coordVert.Num(); i++)
		coordVert[i].Z = coordVert[i].Z * coef;

	return coordVert;
}

void UGenerationFun::ReadVertexArrayFromFile(FString nameFile, TArray<FVector> &VertexArray) {
	TArray<FString> dataArray;
	FString CompleteFilePath = FPaths::ProjectContentDir() + nameFile;
	FFileHelper::LoadANSITextFileToStrings(*CompleteFilePath, NULL, dataArray);
	TArray<FString> currentString;
	FVector currentVector;
	VertexArray.Empty();
	int j;
	for (int i = 0; i < dataArray.Num()-1; i++) {
		dataArray[i].ParseIntoArray(currentString, TEXT(";"), 1);
		//UE_LOG(LogTemp, Warning, TEXT("%d"), currentString.Num());
		j = 0;
		currentVector.X = FCString::Atof(*currentString[j++]);
		currentVector.Y = FCString::Atof(*currentString[j++]);
		currentVector.Z = FCString::Atof(*currentString[j++]);
		VertexArray.Add(currentVector);
	}
}