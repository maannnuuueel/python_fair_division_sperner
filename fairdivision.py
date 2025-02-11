import numpy as np
import matplotlib.pyplot as plt

#implementing the 2-dimensionsal case as described in 
#https://www.cs.cmu.edu/~arielpro/15896s15/docs/paper11a.pdf 

class Triangulation:
    def __init__(self, simplex = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3)/2]])):
        self.points = simplex #2d simplex, triangl
        self.simplices = [[0,1,2]]
        self.labels = [0,1,2] #labeling of points, {0,1,2}, adjacent points have different labels 
        self.colors = []
        
    def bar_coordinates(self,p) -> np.array:
        #euklidian to baycentric coordinates, dim p = 2
        a,b,c = self.points[0], self.points[1], self.points[2]
        T = np.array([[1,1,1], [a[0],b[0],c[0]],[a[1],b[1],c[1]]])
        h = np.array([1,p[0],p[1]])
        out = np.linalg.solve(T,h)
        return out
    
    def eukl_coordinates(self,p) -> np.array:
        #baycentric to euklidian coordinates, dim p = 3
        out = p[0] * self.points[0] + p[1] * self.points[1] + p[2] * self.points[2]
        return out

    def bar_subdivide(self) -> None:
        #barycentric subdivides triangulation 
        self.labels = list(map(lambda x: 0, self.labels)) #everey points in old subdivision gets the same label
        tmp = []
        for s in self.simplices:
            tmp += self.subdivide_simplex(s) #subdive simplex, returns new faces
            self.simplices = [si for si in self.simplices if (s != si)] #remove now subdivided simplex

        self.simplices += tmp

    def subdivide_simplex(self,s) -> list:
        #barycentric subdivides simplex s, returns new faces
        tmp_points_indices=[] #to store indices of added points
        tmp_labels = [1,1,1,2] #labels of added points
        tmp_points = [0.5*(self.points[s[0]] + self.points[s[1]]), 
                      0.5*(self.points[s[1]] + self.points[s[2]]),
                      0.5*(self.points[s[2]] + self.points[s[0]]),
                      (1/3)*(self.points[s[0]] + self.points[s[1]] + self.points[s[2]])]
        
        n=len(self.points) #number of points in triangulation before subdivision

        i,k =0,0 #i counts new points not allready in triangulation, k counts points 
        for point in tmp_points:
            point_not_in_points =1
            
            for j,p in enumerate(self.points): #check if point allready in triangulation
                if (p == point).all(): 
                    tmp_points_indices += [j]
                    point_not_in_points =0
            if point_not_in_points:
                tmp_points_indices += [i+n]
                i += 1
                self.points = np.append(self.points, [point], axis=0)
                self.labels += [tmp_labels[k]] #added points get labeled right away
            k+=1

        return [[s[0],tmp_points_indices[0],tmp_points_indices[3]],[s[0],tmp_points_indices[2],tmp_points_indices[3]],
               [s[1],tmp_points_indices[0],tmp_points_indices[3]],[s[1],tmp_points_indices[1],tmp_points_indices[3]],
               [s[2],tmp_points_indices[1],tmp_points_indices[3]],[s[2],tmp_points_indices[2],tmp_points_indices[3]]]
       
    def plot(self, plot_now = 1, with_color = 0) -> None:
        #plot triangulation
        plt.triplot(self.points[:, 0], self.points[:, 1], self.simplices, color="gray", alpha=0.5)

        #display labels, and colors at points
        markers = ['o','x','v'] #labels
        if with_color:
            colors = ["red", "blue", "green"] #color
            for i, point in enumerate(self.points):
                plt.scatter(point[0], point[1], color=colors[self.colors[i]], marker=markers[self.labels[i]], s=10, label=('room '+ str(self.colors[i]) if i in [0,1,2] else "_"))
        else:
            for i, point in enumerate(self.points):
                plt.scatter(point[0], point[1], color='purple', marker=markers[self.labels[i]], s=10)

        if plot_now:
            plt.title("barycentric triangulation with labels")
            plt.show()  


class RentalHarmony:
    #3 rooms get matched to 3 mates und rent gets divided
    def __init__(self):
        self.tri = Triangulation()
        self.epsilon = 0.0000001

        #not in paper: just to automate to querying
        #self.ideal_shares = [
        #     np.array([1/3, 1/3, 1/3]),  # roommate 0 prefers this
        #     np.array([0.4, 0.3, 0.3]),  # roommate 1 prefers this
        #     np.array([0.35, 0.3, 0.35])   # roommate 2 prefers this
        # ]

    def color_points(self) -> None:
        """Assign colors to each point based on roommate preferences."""
        for i, point in enumerate(self.tri.points):
            self.tri.colors += [self.coloring_fun(point,i)]

    def coloring_fun(self,point,i) -> int:
        rent_division = self.tri.bar_coordinates(point) #barycentric coordinates of point
        
        """if on the outer border color for free room so that coloring is a sperner labeling"""
        zeros = []
        for j,x in enumerate(rent_division):
            if x < self.epsilon/2: zeros += [j]
        
        if len(zeros) == 2: #point is vertex of outer triangle, self.tri.points[0:4]
            return (i + 1) % 3 #just taking the next free one in a cyclic manner so that corners use all the colors
        elif len(zeros) == 1: #point lies on outer edge
            return zeros[0]
        else:
            """for inner points either ask or use ideals"""
            print('rommmate' + str(self.tri.labels[i]) + ': which of the rooms do you prefer in the situation ' + str(rent_division))
            user_input = int(input('choose number 0,1,2:'))
            while user_input not in [0,1,2]: 
                user_input = int(input('you did not enter a number out of 0,1,2, try again:'))
            return user_input
        
            #distances = [rent_division - self.ideal_shares[self.tri.labels[i]]]
            #return int(np.argmin(distances))

    def find_rainbow(self) -> list:
        """find a simplex where every point is colored differently"""
        rainbows=[]
        for simplex in self.tri.simplices:
            simplex_labels = {self.tri.colors[i] for i in simplex}
            if len(simplex_labels) == 3: #number of mates
                rainbows += [simplex]
        return rainbows

    def print_solution(self) -> None:
        #print solution to fair division problem
        rainbows = self.find_rainbow()
        print('Envy-Free Divisions: ')
        for rainbow in rainbows:
            for j in rainbow:
                print(self.tri.bar_coordinates(self.tri.points[j]), 'mate ' + str(self.tri.labels[j]) + ' prefered ' + str(self.tri.colors[j]))
        
    def plot(self) -> None:
        "plot solutions"
        self.tri.plot(plot_now=0, with_color=1)

        #plot self.ideal_shares
        #markers = ['o','x','v'] #labels
        #for i, ideal in enumerate(self.ideal_shares):
        #        tmp = self.tri.eukl_coordinates(ideal)
        #        plt.scatter(tmp[0], tmp[1], marker=markers[i], color="black", s=10, label='ideal roommate ' + str(i))

        # Highlight the rainbow simplex
        rainbows = self.find_rainbow()
        for k,rainbow in enumerate(rainbows):
            rainbow_points = self.tri.points[rainbow]
            plt.fill(rainbow_points[:, 0], rainbow_points[:, 1], "yellow", alpha=0.6, label=('Fair Division' if k==0 else '_'))

        plt.legend()
        plt.title("Rent Division")
        plt.show()


"usage:"

#tri = Triangulation()
#tri.bar_subdivide()
# tri.bar_subdivide()
# tri.bar_subdivide()
# #tri.bar_subdivide()
#tri.plot()

solver = RentalHarmony()
solver.tri.bar_subdivide()
#solver.tri.bar_subdivide()
#solver.tri.bar_subdivide()
#solver.tri.bar_subdivide()
solver.tri.plot()
solver.color_points()
solver.print_solution()
solver.plot()