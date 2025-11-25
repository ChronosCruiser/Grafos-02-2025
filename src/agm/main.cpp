#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <random>
#include <stdexcept>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int DEBUG_OPTION = 0;

class image_load_error
{
};

typedef struct
{
    float peso;
    int vertice1, vertice2;
} aresta;

typedef struct
{
    int rank, peso, tamanho;
} componente;

bool operator<(const aresta &a, const aresta &b)
{
    return a.peso < b.peso;
}

class grafo
{
public:
    grafo(int elementos);
    ~grafo();
    int achar(int vertice);
    void juntar(int vertice1, int vertice2);
    int tamanho(int x);
    int numero_conjuntos();

private:
    int numeros_conjuntos;
    componente *comp;
};

grafo::grafo(int elementos)
{
    comp = new componente[elementos];
    numeros_conjuntos = elementos;

    for (int x = 0; x < elementos; x++)
    {
        comp[x].rank = 0;
        comp[x].tamanho = 1;
        comp[x].peso = x;
    }
}

int grafo::achar(int vertice_index)
{
    int vertice_index_aux = vertice_index;
    while (vertice_index_aux != comp[vertice_index_aux].peso)
    {
        vertice_index_aux = comp[vertice_index_aux].peso;
    }
    comp[vertice_index].peso = vertice_index_aux;
    return vertice_index_aux;
}

void grafo::juntar(int vertice1_index, int vertice2_index)
{
    if (comp[vertice1_index].rank > comp[vertice2_index].rank)
    {
        comp[vertice2_index].peso = vertice1_index;
        comp[vertice1_index].tamanho += comp[vertice2_index].tamanho;
    }
    else
    {
        comp[vertice1_index].peso = vertice2_index;
        comp[vertice2_index].tamanho += comp[vertice1_index].tamanho;
        if (comp[vertice1_index].rank == comp[vertice2_index].rank)
        {
            comp[vertice2_index].rank++;
        }
    }

    numeros_conjuntos--;
}

int grafo::tamanho(int vertice)
{
    return comp[vertice].tamanho;
}

int grafo::numero_conjuntos()
{
    return numeros_conjuntos;
}

struct Pixel
{
    unsigned char red, green, blue;
};

/*
    filename    = nome da imagem para segmentar
    width       = Largura da imagem
    height      = Altura da imagem
    channels    = Numero do canais de cores
*/
std::vector<Pixel> readImage(char *filename, int &width, int &height, int &channels)
{
    unsigned char *data = stbi_load(filename, &width, &height, &channels, 3);
    if (!data)
    {
        throw image_load_error();
    }

    std::vector<Pixel> pixArray(width * height);
    for (int i = 0; i < width * height; i++)
    {
        pixArray[i] = {data[3 * i], data[3 * i + 1], data[3 * i + 2]};
    }

    stbi_image_free(data);
    return pixArray;
}

void debug_print(const char *string)
{
    if (DEBUG_OPTION)
    {
        std::cout << string << std::endl;
    }
}

std::vector<float> fazer_mascara(float sigma)
{
    sigma = std::max(sigma, 0.01F);
    int tamanho = (int)ceil(sigma * 4.0) + 1;
    std::vector<float> mask(tamanho);

    debug_print("[fazer_mascara]: iterando sobre mascara");

    for (int i = 0; i < tamanho; i++)
    {
        mask[i] = exp(-0.5 * (i * i) / (sigma * sigma));
    }

    return mask;
}

/*
    Normalizar a mascara significa aproximar todos os valores entre 0 e 1
    neste caso, usamos uma mascara, normalizando todos od valores da mascara
*/
void normalizar_mascara(std::vector<float> &mascara)
{
    float soma = 0.0f;
    for (size_t i = 0; i < mascara.size(); i++)
    {
        soma += mascara[i];
    }
    soma = 2 * soma - mascara[0];

    for (size_t i = 0; i < mascara.size(); i++)
    {
        mascara[i] /= soma;
    }
}

static void convolucionar_horizontal(const std::vector<float> &entrada, std::vector<float> &resultado,
                                     const std::vector<float> &mascara, int largura, int altura)
{
    int tamanho = mascara.size();

    debug_print("[convolucionar_horizontal]: Convolucionando RGB...");

    for (int y = 0; y < altura; y++)
    {
        for (int x = 0; x < largura; x++)
        {
            float sum = mascara[0] * entrada[y * largura + x];
            for (int i = 1; i < tamanho; i++)
            {
                int pixel_esquerda = std::max(x - i, 0);
                int pixel_direta = std::min(x + i, largura - 1);
                sum += mascara[i] * (entrada[y * largura + pixel_esquerda] + entrada[y * largura + pixel_direta]);
            }
            resultado[y * largura + x] = sum;
        }
    }

    debug_print("[convolucionar_horizontal]: Convolucionando finalizado...");
}

static void convolucionar_vertical(const std::vector<float> &entrada, std::vector<float> &resultado,
                                   const std::vector<float> &mascara, int largura, int altura)
{
    int tamanho = mascara.size();

    debug_print("[convolucionar_vertical]: Convolucionando RGB...");

    for (int x = 0; x < largura; x++)
    {
        for (int y = 0; y < altura; y++)
        {
            float sum = mascara[0] * entrada[y * largura + x];
            for (int i = 1; i < tamanho; i++)
            {
                int top = std::max(y - i, 0);
                int bottom = std::min(y + i, altura - 1);
                sum += mascara[i] * (entrada[top * largura + x] + entrada[bottom * largura + x]);
            }
            resultado[y * largura + x] = sum;
        }
    }

    debug_print("[convolucionar_vertical]: Convolucionando finalizado...");
}


std::vector<float> fazer_gauss(const std::vector<float> &entrada, float sigma, int largura, int altura)
{

    debug_print("[fazer_gauss]: fazendo mascara");
    std::vector<float> mascara = fazer_mascara(sigma);
    debug_print("[fazer_gauss]: normalizar mascara");
    normalizar_mascara(mascara);

    debug_print("[fazer_gauss]: Criar Vetores auxiliares...");
    std::vector<float> temp(largura * altura);
    std::vector<float> saida(largura * altura);

    debug_print("[fazer_gauss]: Convolucionando Vetores...");
    convolucionar_horizontal(entrada, temp, mascara, largura, altura);
    convolucionar_vertical(temp, saida, mascara, largura, altura);

    return saida;
}

/* 
    A dissimilaridade entre pixels sera o metodo principal de verificar se alquer pixel faz parte de um componente
    sendo que a diferenca entre pixels do mesmo componente deve ser baixa e a dissimilaridade entre pixels entre
    componentes diferente deve ser alta.

    Ou seja, pixels do mesmo componente deve ter arestas de peso menor, e pixels de componentes diferentes devem ter
    pesos maior
*/
float dissimilidade(const std::vector<float> &red, const std::vector<float> &green,
                    const std::vector<float> &blue, int largura, int x1, int y1, int x2, int y2)
{
    int idx1 = y1 * largura + x1;
    int idx2 = y2 * largura + x2;
    return sqrt(
        pow(red[idx1] - red[idx2], 2) +
        pow(green[idx1] - green[idx2], 2) +
        pow(blue[idx1] - blue[idx2], 2));
}

grafo *segmentar_grafo(int vertices, aresta *arestas, int numero_arestas, float k)
{
    debug_print("[segmentar_grafo]: Ordenando as arestas...");
    // Como o algoritimo esta usando o Kruskal, deve-se ordernar todas as arestas do menor peso, ate o maior peso
    std::sort(arestas, arestas + numero_arestas);
    grafo *g = new grafo(vertices);

    debug_print("[segmentar_grafo]: Calculando threshold...");
    std::vector<float> threshold(vertices, k);

    debug_print("[segmentar_grafo]: Jutando as arestas...");
    // Para juntar as arestas sera iterado sobre o array prinpical juntando as arestas que tem o menor peso entre elas.
    // Se for achado um peso que esta dentro do nivel minimo sera juntado para um componente apenas.
    for (int i = 0; i < numero_arestas; i++)
    {
        aresta &a = arestas[i];
        int vertice1 = g->achar(a.vertice1);
        int vertice2 = g->achar(a.vertice2);

        if (vertice1 != vertice2 && a.peso <= threshold[vertice1] && a.peso <= threshold[vertice2])
        {
            g->juntar(vertice1, vertice2);
            // Parte de o union find procura um pai comun entre os vertices e junta o novo vertice a ele
            int novo_pai = g->achar(vertice1);
            threshold[novo_pai] = a.peso + k / g->tamanho(novo_pai);
        }
    }
    
    // Retorna a floresta 
    return g;
}

std::vector<Pixel> segmentar_imagem(std::vector<Pixel> imagem, float sigma, float k,
                                    int *num_components, int altura, int largura)
{
    debug_print("[segmentar_imagem]:inicializando os valores de cores....");

    std::vector<float> r(largura * altura);
    std::vector<float> g(largura * altura);
    std::vector<float> b(largura * altura);

    debug_print("[segmentar_imagem]:Copiando os valores da imagem....");

    for (int i = 0; i < largura * altura; i++)
    {
        r[i] = imagem[i].red;
        g[i] = imagem[i].green;
        b[i] = imagem[i].blue;
    }

    debug_print("[segmentar_imagem]:normalizando cores....");
    /*
        Em imagens colorias precisamos diminuir a intencidade de cada coloracao de pixels.
        Por isso usamos o metodo de Gauss para fazer com que os valores dos pixels possam ser
        pesados com mais acuracia, ja que ele pode diminuir os artefatos em uma imagem.
    */
    std::vector<float> smooth_r = fazer_gauss(r, sigma, largura, altura);
    std::vector<float> smooth_g = fazer_gauss(g, sigma, largura, altura);
    std::vector<float> smooth_b = fazer_gauss(b, sigma, largura, altura);

    debug_print("[segmentar_imagem]: criando arestas....");

    /*
        Aqui vamos criar um vetor com todas as arestas possiveis para o nosso grafo,
        ja que cada vertice (pixel) faz vizinha com outro 4 pixeis.
        Sera criado um vetor com todos os possiveis arestas nao direcionadas para cada pixel

        Obs: Pixels perto das bordas nao iram ter 4 vizinhos, podendo apenas 3 ou 2 pixels vizinhos
    */
    std::vector<aresta> arestas;
    arestas.reserve(largura * altura * 4);

    debug_print("[segmentar_imagem]: aplicando peso as arestas arestas....");

    for (int y = 0; y < altura; y++)
    {
        for (int x = 0; x < largura; x++)
        {
            // Fazendo a checagem para ver se o pixel nao vai alem da largura total da imagem
            if (x < largura - 1)
            {
                aresta a;
                a.vertice1 = y * largura + x;
                a.vertice2 = y * largura + (x + 1);
                a.peso = dissimilidade(smooth_r, smooth_g, smooth_b, largura, x, y, x + 1, y);
                arestas.push_back(a);
            }
            // Fazendo a checagem para ver se o pixel nao vai alem da altura total da imagem

            if (y < altura - 1)
            {
                aresta a;
                a.vertice1 = y * largura + x;
                a.vertice2 = (y + 1) * largura + x;
                a.peso = dissimilidade(smooth_r, smooth_g, smooth_b, largura, x, y, x, y + 1);
                arestas.push_back(a);
            }
        }
    }

    
    grafo *graph = segmentar_grafo(largura * altura, arestas.data(), arestas.size(), k);
    *num_components = graph->numero_conjuntos();

    /* 
        Metodos para escolher e criar cores aleatorias com base no tamanho de componentes achados
    */
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<unsigned char> dist(0, 255);
    std::unordered_map<int, Pixel> cores;

    for (int y = 0; y < altura; y++)
    {
        for (int x = 0; x < largura; x++)
        {
            int idx = y * largura + x;
            int pai = graph->achar(idx);
            if (cores.find(pai) == cores.end())
            {
                //intera sobre as cores criando cores aleatorias
                cores[pai] = {dist(rng), dist(rng), dist(rng)};
            }
        }
    }

    /* 
     Cria uma nova imagem com os pixels e segmentacoes
    */
    std::vector<Pixel> output(largura * altura);
    for (int y = 0; y < altura; y++)
    {
        for (int x = 0; x < largura; x++)
        {
            int idx = y * largura + x;
            output[idx] = cores[graph->achar(idx)];
        }
    }

    // delete graph;
    return output;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cerr << "Uso: " << argv[0] << " <caminho da imagem> <caminho da saida>\n";
        return 1;
    }

    int width, height, channels;
    std::vector<Pixel> img = readImage(argv[1], width, height, channels);

    /*
        valores padroes para ser usado no algoritimo
        sigma = taxa de nivelamento de cores de pixel antes de fazer o grafo
        k = valor para a funcao de teto maximo
        components = Numero de componentes total da imagem
    */
    float sigma = 0.5;
    float k = 800;
    int components = 0;

    std::vector<Pixel> segmented = segmentar_imagem(img, sigma, k, &components, height, width);
    stbi_write_png(argv[2], width, height, 3, segmented.data(), width * 3);

    std::cout << "Segmentacao completa. Achou " << components << " componentes.\n";

    return 0;
}