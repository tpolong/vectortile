package main

import (
	"encoding/csv"
	"errors"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"net/http"
	"regexp"
	"strconv"
	"strings"
	"time"
	//"math/rand"
	"github.com/golang/protobuf/proto"
	"github.com/vicapow/go-vtile-example/gen/third-party/vector-tile-spec/1.0.1"
)

func cmdEnc(id uint32, count uint32) uint32 {
	return (id & 0x7) | (count << 3)
}

func moveTo(count uint32) uint32 {
	return cmdEnc(1, count)
}

func lineTo(count uint32) uint32 {
	return cmdEnc(2, count)
}

func closePath(count uint32) uint32 {
	return cmdEnc(7, count)
}

func paramEnc(value int32) int32 {
	return (value << 1) ^ (value >> 31)
}
func tile2lon( x int,  z int)(a float64) {
	return float64(x) /math.Pow(2, float64(z)) * 360.0 - 180;
 }

 func tile2lat( y int,  z int)(a float64) {
   n := math.Pi - (2.0 * math.Pi * float64(y)) / math.Pow(2, float64(z));
   return math.Atan(math.Sinh(n))*180/math.Pi;
 }
 
func FloatToString(input_num float64) string {
    // to convert a float number to a string
    return strconv.FormatFloat(input_num, 'f', 6, 64)
}
func createTileWithPoints(xyz TileID, pnt []Point,brush *KDBush) ([]byte, error) {
	tile := &vector_tile.Tile{}
	var layerVersion = vector_tile.Default_Tile_Layer_Version
	layerName := "points"
	featureType := vector_tile.Tile_POINT
	var extent = vector_tile.Default_Tile_Layer_Extent
	// Put a point in the center of the tile.
	var geometry []uint32
	var filtered [][2]float64
	ymax :=tile2lat(int(xyz.y), int(xyz.z));
	ymin := tile2lat(int(xyz.y+1), int(xyz.z));
	xmin := tile2lon(int(xyz.x), int(xyz.z));
	xmax :=tile2lon(int(xyz.x+1), int(xyz.z));
	// fmt.Println("ymax: ", ymax)
	// fmt.Println("ymin: ",ymin)
	// fmt.Println("xmin : ",xmin )
	// fmt.Println("xmax : ",xmax )
	// for _, point := range points {
	// 	x, y := lngLatToTileXY(point, xyz)
	// 	if x >= 0 && x < 1 && y >= 0 && y < 1 {
	// 		filtered = append(filtered, [2]float64{x, y})
	// 	}
	// }
	
	result :=  brush.Range(xmin ,ymin, xmax, ymax)
	for _, i := range result{
		a,b:=pnt[i].Coordinates()
		x, y := lngLatToTileXY(LngLat{lng:a,lat:b}, xyz)
		if x >= 0 && x < 1 && y >= 0 && y < 1 {
			filtered = append(filtered, [2]float64{x, y})
		}
	}
	//fmt.Println("result : ", result )
	if len(filtered) > 0 {
		cmd := moveTo(uint32(len(filtered)))
		geometry = append(geometry, cmd)
		var pX int32
		var pY int32
		for _, point := range filtered {
			deltaX := int32(float64(extent)*point[0]+0.5) - pX
			deltaY := int32(float64(extent)*point[1]+0.5) - pY
			geometry = append(geometry, uint32(paramEnc(deltaX)))
			geometry = append(geometry, uint32(paramEnc(deltaY)))
			pX = pX + deltaX
			pY = pY + deltaY
		}
	} else {
		// Return an empty tile if we have no points
		return nil, nil
	}
	tile.Layers = []*vector_tile.Tile_Layer{
		&vector_tile.Tile_Layer{
			Version: &layerVersion,
			Name:    &layerName,
			Extent:  &extent,
			Features: []*vector_tile.Tile_Feature{
				&vector_tile.Tile_Feature{
					Tags:     []uint32{},
					Type:     &featureType,
					Geometry: geometry,
				},
			},
		},
	}
	return proto.Marshal(tile)
}

func xyzToLngLat(tileX float64, tileY float64, tileZ float64) (float64, float64) {
	totalTilesX := math.Pow(2, tileZ)
	totalTilesY := math.Pow(2, tileZ)
	x := float64(tileX) / float64(totalTilesX)
	y := float64(tileY) / float64(totalTilesY)
	// lambda can go from [-pi/2, pi/2]
	lambda := x*math.Pi*2 - math.Pi
	// phi can go from [-1.4844, 1.4844]
	phi := 2*math.Atan(math.Exp((2*y-1)*math.Pi)) - (math.Pi / 2)
	lng := lambda * 180 / math.Pi
	lat := (math.Pi - phi) * 180 / math.Pi
	return lng, lat
}

func lngLatToTileXY(ll LngLat, tile TileID) (float64, float64) {
	totalTilesX := math.Pow(2, float64(tile.z))
	totalTilesY := math.Pow(2, float64(tile.z))
	lambda := (ll.lng + 180) / 180 * math.Pi
	// phi: [-pi/2, pi/2]
	phi := ll.lat / 180 * math.Pi
	tileX := lambda / (2 * math.Pi) * totalTilesX
	// [-1.4844, 1.4844] -> [1, 0]  * totalTilesY
	tileY := (math.Log(math.Tan(math.Pi/4-phi/2))/math.Pi/2 + 0.5) * totalTilesY
	return tileX - float64(tile.x), tileY - float64(tile.y)
}

// Takes a string of the form `<z>/<x>/<y>` (for example, 1/2/3) and returns
// the individual uint32 values for x, y, and z if there was no error.
// Otherwise, err is set to a non `nil` value and x, y, z are set to 0.
func tilePathToXYZ(path string) (TileID, error) {
	xyzReg := regexp.MustCompile("(?P<z>[0-9]+)/(?P<x>[0-9]+)/(?P<y>[0-9]+)")
	matches := xyzReg.FindStringSubmatch(path)
	if len(matches) == 0 {
		return TileID{}, errors.New("Unable to parse path as tile")
	}
	x, err := strconv.ParseUint(matches[2], 10, 32)
	if err != nil {
		return TileID{}, err
	}
	y, err := strconv.ParseUint(matches[3], 10, 32)
	if err != nil {
		return TileID{}, err
	}
	z, err := strconv.ParseUint(matches[1], 10, 32)
	if err != nil {
		return TileID{}, err
	}
	return TileID{x: uint32(x), y: uint32(y), z: uint32(z)}, nil
}

// A LngLat is a struct that holds a longitude and latitude vale.
type LngLat struct {
	lng float64
	lat float64
}

// TileID represents the id of the tile.
type TileID struct {
	x uint32
	y uint32
	z uint32
}

// Tree a struct holder for tree information.
// type Tree struct {
// 	lng     float64
// 	lat     float64
// 	species string
// }
type Tree struct {
	lng     float64
	lat     float64
	//species string
}
func loadTrees() []Point{
	content, err := ioutil.ReadFile("point100w.csv")
	if err != nil {
		log.Fatal(err)
	}
	r := csv.NewReader(strings.NewReader(string(content[:])))
	records, err := r.ReadAll()
	if err != nil {
		log.Fatal(err)
	}
	// TreeID,qLegalStatus,qSpecies,qAddress,SiteOrder,qSiteInfo,PlantType,qCaretaker,
	// qCareAssistant,PlantDate,DBH,PlotSize,PermitNotes,XCoord,YCoord,Latitude,Longitude,Location
	//var pnt []Point
	pnt:= make([]Point, len(records[0:]))
	//fmt.Println("数据: ", )
	var i int64=0
	for _, record := range records[0:]{

		lat, _ := strconv.ParseFloat(record[0], 64)
		lng, _ := strconv.ParseFloat(record[1], 64)
		//trees = append(trees, Tree{lng: lng, lat: lat})
		pnt[i]=&SimplePoint{X:lng,Y:lat}
		i=i+1
	}
	return pnt
}
// Interface, that should be implemented by indexing structure
// It's just simply returns points coordinates
// Called once, only when index created, so you could calc values on the fly for this interface
type Point interface {
	Coordinates() (X, Y float64)
}

// Minimal struct, that implements Point interface
type SimplePoint struct {
	X, Y float64
}

// SimplePoint's  implementation of Point interface
func (sp *SimplePoint) Coordinates() (float64, float64) {
	return sp.X, sp.Y
}

// A very fast static spatial index for 2D points based on a flat KD-tree.
// Points only, no rectangles
// static (no add, remove items)
// 2 dimensional
// indexing 16-40 times faster then  rtreego(https://github.com/dhconnelly/rtreego) (TODO: benchmark)
type KDBush struct {
	NodeSize int
	Points   []Point

	idxs   []int     //array of indexes
	coords []float64 //array of coordinates
}

// Create new index from points
// Structure don't copy points itself, copy only coordinates
// Returns pointer to new KDBush index object, all data in it already indexed
// Input:
// points - slice of objects, that implements Point interface
// nodeSize  - size of the KD-tree node, 64 by default. Higher means faster indexing but slower search, and vise versa.
func NewBush(points []Point, nodeSize int) *KDBush {
	b := KDBush{}
	b.buildIndex(points, nodeSize)
	return &b
}

// Finds all items within the given bounding box and returns an array of indices that refer to the items in the original points input slice.
func (bush *KDBush) Range(minX, minY, maxX, maxY float64) []int {
	stack := []int{0, len(bush.idxs) - 1, 0}
	result := []int{}
	var x, y float64

	for len(stack) > 0 {
		axis := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		right := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		left := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		if right-left <= bush.NodeSize {
			for i := left; i <= right; i++ {
				x = bush.coords[2*i]
				y = bush.coords[2*i+1]
				if x >= minX && x <= maxX && y >= minY && y <= maxY {
					result = append(result, bush.idxs[i])
				}
			}
			continue
		}

		m := floor(float64(left+right) / 2.0)

		x = bush.coords[2*m]
		y = bush.coords[2*m+1]

		if x >= minX && x <= maxX && y >= minY && y <= maxY {
			result = append(result, bush.idxs[m])
		}

		nextAxis := (axis + 1) % 2

		if (axis == 0 && minX <= x) || (axis != 0 && minY <= y) {
			stack = append(stack, left)
			stack = append(stack, m-1)
			stack = append(stack, nextAxis)
		}

		if (axis == 0 && maxX >= x) || (axis != 0 && maxY >= y) {
			stack = append(stack, m+1)
			stack = append(stack, right)
			stack = append(stack, nextAxis)
		}

	}
	return result
}

// Finds all items within a given radius from the query point and returns an array of indices.
func (bush *KDBush) Within(point Point, radius float64) []int {
	stack := []int{0, len(bush.idxs) - 1, 0}
	result := []int{}
	r2 := radius * radius
	qx, qy := point.Coordinates()

	for len(stack) > 0 {
		axis := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		right := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		left := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		if right-left <= bush.NodeSize {
			for i := left; i <= right; i++ {
				dst := sqrtDist(bush.coords[2*i], bush.coords[2*i+1], qx, qy)
				if dst <= r2 {
					result = append(result, bush.idxs[i])
				}
			}
			continue
		}

		m := floor(float64(left+right) / 2.0)
		x := bush.coords[2*m]
		y := bush.coords[2*m+1]

		if sqrtDist(x, y, qx, qy) <= r2 {
			result = append(result, bush.idxs[m])
		}

		nextAxis := (axis + 1) % 2

		if (axis == 0 && (qx-radius <= x)) || (axis != 0 && (qy-radius <= y)) {
			stack = append(stack, left)
			stack = append(stack, m-1)
			stack = append(stack, nextAxis)
		}

		if (axis == 0 && (qx+radius >= x)) || (axis != 0 && (qy+radius >= y)) {
			stack = append(stack, m+1)
			stack = append(stack, right)
			stack = append(stack, nextAxis)
		}
	}

	return result
}

///// private method to sort the data

////////////////////////////////////////////////////////////////
/// Sorting stuff
////////////////////////////////////////////////////////////////

func (bush *KDBush) buildIndex(points []Point, nodeSize int) {
	bush.NodeSize = nodeSize
	bush.Points = points

	bush.idxs = make([]int, len(points))
	bush.coords = make([]float64, 2*len(points))

	for i, v := range points {
		bush.idxs[i] = i
		x, y := v.Coordinates()
		bush.coords[i*2] = x
		bush.coords[i*2+1] = y
	}

	sort(bush.idxs, bush.coords, bush.NodeSize, 0, len(bush.idxs)-1, 0)
}

func sort(idxs []int, coords []float64, nodeSize int, left, right, depth int) {
	if (right - left) <= nodeSize {
		return
	}

	m := floor(float64(left+right) / 2.0)

	sselect(idxs, coords, m, left, right, depth%2)

	sort(idxs, coords, nodeSize, left, m-1, depth+1)
	sort(idxs, coords, nodeSize, m+1, right, depth+1)

}

func sselect(idxs []int, coords []float64, k, left, right, inc int) {
	//whatever you want
	for right > left {
		if (right - left) > 600 {
			n := right - left + 1
			m := k - left + 1
			z := math.Log(float64(n))
			s := 0.5 * math.Exp(2.0*z/3.0)
			sds := 1.0
			if float64(m)-float64(n)/2.0 < 0 {
				sds = -1.0
			}
			n_s := float64(n) - s
			sd := 0.5 * math.Sqrt(z*s*n_s/float64(n)) * sds
			newLeft := iMax(left, floor(float64(k)-float64(m)*s/float64(n)+sd))
			newRight := iMin(right, floor(float64(k)+float64(n-m)*s/float64(n)+sd))
			sselect(idxs, coords, k, newLeft, newRight, inc)
		}

		t := coords[2*k+inc]
		i := left
		j := right

		swapItem(idxs, coords, left, k)
		if coords[2*right+inc] > t {
			swapItem(idxs, coords, left, right)
		}

		for i < j {
			swapItem(idxs, coords, i, j)
			i += 1
			j -= 1
			for coords[2*i+inc] < t {
				i += 1
			}
			for coords[2*j+inc] > t {
				j -= 1
			}
		}

		if coords[2*left+inc] == t {
			swapItem(idxs, coords, left, j)
		} else {
			j += 1
			swapItem(idxs, coords, j, right)
		}

		if j <= k {
			left = j + 1
		}
		if k <= j {
			right = j - 1
		}
	}
}

func swapItem(idxs []int, coords []float64, i, j int) {
	swapi(idxs, i, j)
	swapf(coords, 2*i, 2*j)
	swapf(coords, 2*i+1, 2*j+1)
}

func swapf(a []float64, i, j int) {
	t := a[i]
	a[i] = a[j]
	a[j] = t
}

func swapi(a []int, i, j int) {
	t := a[i]
	a[i] = a[j]
	a[j] = t
}

func iMax(a, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}

func iMin(a, b int) int {
	if a < b {
		return a
	} else {
		return b
	}
}

func floor(in float64) int {
	out := math.Floor(in)
	return int(out)
}

func sqrtDist(ax, ay, bx, by float64) float64 {
	dx := ax - bx
	dy := ay - by
	return dx*dx + dy*dy
}

func main() {
	t1 := time.Now() 
	pnt := loadTrees()
	//points := make([]LngLat, len(trees))
	// pnt := make([]Point, len(trees))
	// for i, tree := range trees {
	// 	//points[i] = LngLat{lng: tree.lng, lat: tree.lat}
	// 	pnt[i]=&SimplePoint{X:tree.lng,Y:tree.lat}
	// }
	// points = points[0:1000]
	bush := NewBush(pnt, 10)
	//fmt.Println("number of points", len(points))
	mux := http.NewServeMux()
	elapsed := time.Since(t1)
	fmt.Println("计时: ", elapsed)
	// Handle requests for urls of the form `/tiles/{z}/{x}/{y}` and returns
	// the vector tile for the even tile x, y, and z coordinates.
	tileBase := "/tiles/"
	mux.HandleFunc(tileBase, func(w http.ResponseWriter, r *http.Request) {
		t2 := time.Now() 
		log.Printf("url: %s", r.URL.Path)
		tilePart := r.URL.Path[len(tileBase):]
		xyz, err := tilePathToXYZ(tilePart)
		if err != nil {
			http.Error(w, "Invalid tile url", 400)
			return
		}
		data, err := createTileWithPoints(xyz, pnt,bush)
		if err != nil {
			log.Fatal("error generating tile", err)
		}
		// All this APi to be requests from other domains.
		w.Header().Set("Content-Type", "application/x-protobuf")
		w.Header().Set("Access-Control-Allow-Origin", "*")
		w.Header().Set("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
		w.Write(data)
		elapsed2 := time.Since(t2)
	
	fmt.Println("耗时: ", elapsed2)
	})
	log.Fatal(http.ListenAndServe(":8081", mux))
}
