from flask import Flask, request, jsonify
import tensorflow as tf
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split

app = Flask(__name__)

model1 = tf.keras.models.load_model("E:/New_P/NEW/flask_h5/BEM-GRUcnnf.h5")
model2 = tf.keras.models.load_model("E:/New_P/NEW/flask_h5/S175_BEM-GRUcnnf.h5")
model3 = tf.keras.models.load_model("E:/New_P/NEW/flask_h5/s60_BEM-GRUcnnf .h5")
model4 = tf.keras.models.load_model("E:/New_P/NEW/flask_h5/wigleyBEM-GRUcnnf.h5")

@app.route

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()  
    input_data = np.array(data['input']) 
    input_data1 = input_data.reshape(1, 7)
    input_data_reshaped = input_data1.reshape(-1, 1, 7)
    prediction = model1.predict(input_data_reshaped)
    prediction1 = prediction
    prediction = model2.predict(input_data_reshaped)
    prediction2 = prediction
    prediction = model3.predict(input_data_reshaped)
    prediction3 = prediction
    prediction = model4.predict(input_data_reshaped)
    prediction4 = prediction

    predictions_matrix = np.vstack((prediction1, prediction2, prediction3, prediction4))
    return jsonify({'predictions': predictions_matrix.tolist()})
if __name__ == '__main__':
    app.run(debug=True)

