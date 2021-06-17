#include "widget.h"
#include "ui_widget.h"
#include "mathematics.cpp"

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
    ui->setupUi(this);
    scene = new QGraphicsScene(ui->graphicsView);
    QPen pen(Qt::black);
    scene->addLine((0)*100+5,-(0)*100+5,5,-1000,pen);
    scene->addLine(5,-1000,15,-980,pen);
    scene->addLine(5,-1000,-5,-980,pen);
    scene->addLine((0)*100+5,-(0)*100+5,1000,5,pen);
    scene->addLine(980,-5,1000,5,pen);
    scene->addLine(980,15,1000,5,pen);
    text = new  QGraphicsTextItem ;
}

Widget::~Widget()
{
    delete ui;
}



void Widget::on_pushButton_clicked()
{
    QPen pen(Qt::white);
    pen.setWidth(3);
    QPen pen1(Qt::blue);
    QBrush br(QColor(220,0,0));
    QBrush br1(QColor(255,0,0));

    double W_graph;
    double T_graph;
    table(W_graph, T_graph);
    //if (x->back() ==0)
    //{
    //    scene->addLine((T_graph)*100+5,-(W_graph)*100+5,x->back()+5,y->back()+5,pen1);
    //}
    scene->addEllipse((T_graph)*100,-(W_graph)*100,10,10,pen,br);
    //scene->addLine((T_graph)*100+5,-(W_graph)*100+5,x->last()+5,y->last()+5,pen1);
    text->setX((T_graph) *100 -10);
    text->setY((-W_graph)*100 +10);
    text->setPlainText("("+QString::number(T_graph)+","+QString::number(W_graph)+")");
    //scene->addText();
    //scene->addText("("+QString::number(T_graph)+","+QString::number(W_graph)+")");
    scene->addItem(text);
    //y->push_back(W_graph);
    //x->push_back(T_graph);
    ui->graphicsView->setScene(scene);
}
