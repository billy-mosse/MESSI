import traceback

from flask import Flask
from flask import request, abort, make_response

import json
app = Flask(__name__)
app.config["TRAP_HTTP_EXCEPTIONS"] = True


@app.route("/process", methods = ['POST'])
def process_pdf():
    net = request.form.get("net").replace('\n', '')
    return json.dumps({'answer': net})


@app.errorhandler(404)
def resource_not_found(e):
    return {"error": "resource not found", "details": str(e)}, 404


@app.errorhandler(Exception)
def global_handler(e):
    traceback.print_exc()
    return {"error": "something bad happened", "details": str(e)}, 500





if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)
