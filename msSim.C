#include <iostream>
#include <TROOT.h>
#include "math.h"
#include <random>

// g++ -o msSim -c msSim.C
/*
def createScatteringRanges(low,high):

    numScatter =  random.randint(low,high)
    begin, end = 0, 1
    for j in range(numScatter):
        msScatterDistance.append(random.uniform(begin,end))
        #begin = msScatterDistance[j]
    #print msScatterDistance
    #print random.gauss(0, 0.1)
    msScatterDistance.sort()
*/
std::vector<float> createScatteringRanges(std::exponential_distribution<> rng, std::mt19937 & rnd_gen , float pathLength){

    std::vector<float> scatter_distances;
    float scatter_distance;
    float total_distance = 0;
    while(total_distance < pathLength*1.5){
        scatter_distance = (float)rng (rnd_gen) * 4;
        
        scatter_distances.push_back(scatter_distance);

        total_distance = total_distance + scatter_distance;
    
    }

    return scatter_distances;
}

std::vector<float> createScatteringAngles(std::vector<float> scatter_vector, float scatter_angle, std::mt19937 & rnd_gen, std::normal_distribution<double> scatter_dist){

    std::vector<float> scattering_angles_vector;

    for(int i =0; i< scattering_angles_vector.size(); i++){

         std::cout << scatter_dist(rnd_gen) << std::endl;

    }

    return scattering_angles_vector;
}


int msSim(){
    float x_pos_detector = 1.0;
    float y_target_pos  = 1.0;
    float pathLength = sqrt(x_pos_detector*x_pos_detector + y_target_pos*y_target_pos);

    int num = 10;
    TH1F TH1F_x_residual = TH1F("TH1F_x_residual", "x residual at ideal position; x; counts", 100, -.5,.5);
    TH1F TH1F_x_residual_corrected = TH1F("TH1F_x_residual_corrected", "x residual at ideal position corrected; x; counts", 100, -.5,.5);
    TH2F TH2F_x_residual_dx_residual = TH2F("TH2F_x_residual_dx_residual", "x vs dx/dy residual", 100, -1,1,100, 1,1);
    TH2F TH2F_x_residual_v_ms = TH2F("TH2F_x_residual_v_ms", "x vs dx/dy residual", 100, -1,1,100, 0,100);


    float scatterLengthScale = x_pos_detector/80;
    float scatterLengthScaleLambda = 1/scatterLengthScale;
    float scatter_angle = .01;

    std::random_device rd; 
    std::mt19937 rnd_gen (rd());

    std::exponential_distribution<> rng (scatterLengthScaleLambda);
    std::normal_distribution<double> scatter_angle_dist(0.0,scatter_angle);

    
    std::vector<float> scatter_vector;

    for(int i=0;i<num;i++){
/*
        scatter_vector = createScatteringRanges(rng,rnd_gen, pathLength);
        createScatteringAngles(scatter_vector, scatter_angle, rnd_gen, scatter_dist);
*/

        std::vector<float> scatter_distances;
        float scatter_distance;
        float total_distance = 0;
        std::vector<float> scattering_angles_vector;
        float x_pos = 0;
        float y_pos = 0;
        float scatter_angle_l;
        while(total_distance < pathLength*1.5){
            scatter_distance = (float)rng (rnd_gen) * 4;

            scatter_angle_l = scatter_angle_dist(rnd_gen);

            std::cout << scatter_angle_l << " " << scatter_distance << " " << x_pos << " "  << y_pos  << std::endl;
    
            total_distance = total_distance + scatter_distance;
            x_pos = scatter_distance*cos(scatter_angle_l) + x_pos;
            y_pos = scatter_distance*sin(scatter_angle_l) + x_pos;
        
        }

        
        if(i*10 % num == 0){
            std::cout <<  (float)i/(float)num << std::endl;
            
        }
    
    }
    return 0;
}

int main(){
    

std::cout <<  "test" << std::endl;

/*
for i in range(num):
    if i%1000 == 0: print (i+ 0.0)/num
    m, b = target, 0.0
    msScatterDistance= []
    line = [1,0]
    createScatteringRanges(25,100)
    #print "new event" , msScatterDistance
    line = scatterTrack(msScatterDistance,m,b, 0.01)
    #print line
    residual_Value = 0;
    if line: 
        #print  "residual" , residual(line[0],line[1]), line[0], line[1]
        residual_Value = residual(line[0],line[1]) - target
    #print "residual" , residual_Value
    TH1F_x_residual.Fill(residual_Value)
    phi = math.atan((residual_Value + target)/1.0) - math.atan(target/1.0)
    correctedX = phi*math.sqrt((1 + math.pow(residual_Value -target/4, 2)))
    dxdyRes = target-line[0]
    correctedX = residual_Value+ .6*dxdyRes
    TH1F_x_residual_corrected.Fill(correctedX)

    #print dxdyRes

    #TH2F_x_residual_dx_residual.Fill(correctedX, dxdyRes)
    TH2F_x_residual_dx_residual.Fill(residual_Value, dxdyRes)
    TH2F_x_residual_v_ms.Fill(residual_Value, len(msScatterDistance))
*/

}


