// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "GenerationFun.generated.h"

/**
 * 
 */
UCLASS()
class STRATEGYGAME_API UGenerationFun : public UBlueprintFunctionLibrary
{
	GENERATED_BODY() public:
			UFUNCTION(BlueprintCallable, Category = "Cutom", meta = (Keywords = "Save data in file"))
				static void SaveDataToFile(const TArray<FString>& Data, FString Name);
			UFUNCTION(BlueprintCallable, Category = "Cutom", meta = (Keywords = "NextString"))
				static FString NextString(FString Data);
			UFUNCTION(BlueprintCallable, Category = "Custom", meta = (Keywords = "Read Vertex Array From File"))
				static void ReadVertexArrayFromFile(FString nameFile, TArray<FVector>& VertexArray);

			UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "ErosionScore"))
				static float erosionScore(const TArray<FVector>& coordVert, TArray<float> &matrixS );
			UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Change Vertex Array Coord"))
				static TArray<FVector> changeVertexArrayCoord(const TArray<FVector>& arrayVertex, float coefXY, float coefZ);
			UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Create Vertex Array Component"))
				static TArray<FVector> CreateVertexArrayComponent(const TArray<FVector>& arrayVertex, int indexComponent, int width);
			UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Create Trees"))
				static TArray<FVector> CreateTrees(const TArray<FVector>& coordVert);
			UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Normalize Coord Vector"))
				static TArray<FVector> NormalizeCoordVector(TArray<FVector> coordVert);


			UFUNCTION(BlueprintCallable, Category = "NavMesh")
				static void UpdateNavMesh_Actor(AActor* Actor)
			{
				AActor& AVal = *Actor;
				FNavigationSystem::UpdateActorAndComponentData(AVal, true);
			}

			UFUNCTION(BlueprintCallable, Category = "NavMesh")
				static void UpdateNavMesh_Component(UActorComponent* ActorComponent)
			{
				UActorComponent& CVal = *ActorComponent;
				FNavigationSystem::UpdateComponentData(CVal);
			}
};
